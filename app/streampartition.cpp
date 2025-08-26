/******************************************************************************
 * streampartition.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj <marcelofaraj@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <memory>
#include <optional>
#include <cmath>

#include "definitions.h"
#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/buffer.h"
#include "graph_io_stream.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "partition/uncoarsening/refinement/mixed_refinement.h"
#include "partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "tools/flat_buffer_writer.h"
#include "mlp_thread_manager.cpp"
#include "batch_id_manager.h"


long getMaxRSS();

void debug_print(std::string msg) {
    std::cout << msg << std::endl;
}

std::string extractBaseFilename(const std::string &fullPath);

int main(int argn, char **argv) {
    /* std::cout << R"(
    ██   ██ ███████ ██ ███████ ████████ ██████  ███████  █████  ███    ███
    ██   ██ ██      ██ ██         ██    ██   ██ ██      ██   ██ ████  ████
    ███████ █████   ██ ███████    ██    ██████  █████   ███████ ██ ████ ██
    ██   ██ ██      ██      ██    ██    ██   ██ ██      ██   ██ ██  ██  ██
    ██   ██ ███████ ██ ███████    ██    ██   ██ ███████ ██   ██ ██      ██


    ███    ██  ██████  ██████  ███████
    ████   ██ ██    ██ ██   ██ ██
    ██ ██  ██ ██    ██ ██   ██ █████
    ██  ██ ██ ██    ██ ██   ██ ██
    ██   ████  ██████  ██████  ███████
    )" << std::endl; */

    timer processing_t, io_t, t_mlp, first_phase_t, second_phase_t, part_single_node_t, buffer_add_node_t, wait_for_mlp_t;
    TIMING_DECLARE(double global_mapping_time) = 0;
    TIMING_DECLARE(double buffer_io_time) = 0;
    TIMING_DECLARE(double io_time) = 0;
    TIMING_DECLARE(double buffer_add_node_time) = 0;
    TIMING_DECLARE(double model_construction_time) = 0;
    TIMING_DECLARE(double first_phase_time) = 0;
    TIMING_DECLARE(double second_phase_time) = 0;
    TIMING_DECLARE(double updating_adj_time) = 0;
    TIMING_DECLARE(double part_single_node_time) = 0;
    TIMING_DECLARE(double mlp_time) = 0;
    TIMING_DECLARE(double wait_for_mlp_finish_time) = 0;

    PartitionConfig partition_config;
    std::string graph_filename;
    uint64_t total_edge_cut = 0;

    quality_metrics qm;
    balance_configuration bc;
    std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> *batch_nodes;

    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;

    int ret_code = parse_parameters(argn, argv,
                                    partition_config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    if (ret_code) {
        return 0;
    }

    std::ofstream ofs;
    ofs.open("/dev/null");
    if (suppress_output) {
        std::cout.rdbuf(ofs.rdbuf());
    }
    srand(partition_config.seed);
    random_functions::setSeed(partition_config.seed);

    partition_config.LogDump(stdout);
    partition_config.graph_filename = graph_filename;
    partition_config.stream_input = true;

    partition_config.batch_manager = new BatchIDManager(partition_config.max_active_batches);
    size_t cur_batch_id = partition_config.batch_manager->acquire_id();
    PartitionID cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);
    std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>* current_batch = nullptr; //new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>();

    int &passes = partition_config.num_streams_passes;
    partition_config.count_misc1 = 0;
    partition_config.count_misc2 = 0;

    Buffer* buffer = new Buffer(partition_config, partition_config.max_pq_size);
    MLPThreadManager mlp_thread_manager;

    for (partition_config.restream_number = 0; partition_config.restream_number < passes; partition_config.restream_number++) {

        partition_config.store_unpartitioned_neighbors = partition_config.ghost_importance > 0 && partition_config.restream_number == 0;

        // ***************************** IO operations ***************************************
        TIMING_START(io_t);
        graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);

        double avg_block_size = static_cast<double>(partition_config.number_of_nodes) / partition_config.k;
        partition_config.max_block_weight = static_cast<int>(std::ceil((1.0 + partition_config.imbalance / 100) * avg_block_size));
        TIMING_ACCUMULATE(io_time, io_t);

        if (partition_config.restream_number == 1) {
            // partition_config.stream_buffer_len = MAX(partition_config.max_pq_size / 10, partition_config.stream_buffer_len);
        }

        TIMING_START(first_phase_t);

        std::vector<std::string> lines(1);
        buffered_input ss2(&lines);
        std::vector<LongNodeID> cur_line;
        cur_line.reserve(1000);

        bool use_mlp = partition_config.stream_buffer_len != 1;

        partition_config.total_nodes_loaded = 0;
        while (partition_config.remaining_stream_nodes) {

            TIMING_START(io_t);
            // Load a line from the stream
            std::getline(*(partition_config.stream_in), lines[0]);
            if (lines[0][0] == '%') { // skip comments in the file
                continue;
            }

            partition_config.remaining_stream_nodes--;
            LongNodeID global_node_id = ++partition_config.total_nodes_loaded;

            assert(global_node_id <= partition_config.number_of_nodes);

            ss2.simple_scan_line_fast(cur_line);
            TIMING_ACCUMULATE(io_time, io_t);

            unsigned degree = cur_line.size();

            if (degree >= partition_config.d_max || degree == 0) { //
                // Partition node directly if degree is too high or 0

                if (partition_config.restream_number == 0) {
                    if (degree == 0) {
                        // Put node into partition with lowest weight
                        PartitionID best_partition = 0;
                        LongNodeID min_weight = partition_config.max_block_weight;

                        // Wait for MLP to finish
                        TIMING_START(wait_for_mlp_t);
                        mlp_thread_manager.wait_completion();
                        TIMING_ACCUMULATE(wait_for_mlp_finish_time, wait_for_mlp_t);

                        for (PartitionID i = 0; i < partition_config.k; i++) {
                            if ((*partition_config.stream_blocks_weight)[i] < min_weight) {
                                min_weight = (*partition_config.stream_blocks_weight)[i];
                                best_partition = i;
                            }
                        }
                        (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;
                        (*partition_config.stream_blocks_weight)[best_partition]++;


                    } else {

                        // Partition node directly
                        partition_config.count_misc1++;

                        // Wait for MLP to finish
                        TIMING_START(wait_for_mlp_t);
                        mlp_thread_manager.wait_completion();
                        TIMING_ACCUMULATE(wait_for_mlp_finish_time, wait_for_mlp_t);

                        TIMING_START(part_single_node_t);
                        partition_single_node(partition_config, global_node_id, cur_line);
                        TIMING_ACCUMULATE(part_single_node_time, part_single_node_t);

                        // Update neighbors
                        buffer->update_neighbours_priority(cur_line);
                    }
                    continue;
                } else if (!partition_config.restream_include_high_degree_nodes ) {
                    // In restreaming if not including high degree nodes, we do not process nodes with degree 0 or degree > d_max
                    cur_line.clear();
                    TIMING_ACCUMULATE(first_phase_time, first_phase_t);
                    continue; // Skip empty or too large nodes
                }

            }

            // debug_print(std::to_string((getMaxRSS())));
            if (partition_config.restream_number == 0) {
                // Check if new node has a higher buffer score than max score in buffer
                TIMING_START(buffer_add_node_t);
                bool added_to_buffer = buffer->addNode(global_node_id, cur_line);
                TIMING_ACCUMULATE(buffer_add_node_time, buffer_add_node_t);

                if (!added_to_buffer) {
                    // Wait for MLP to finish
                    TIMING_START(wait_for_mlp_t);
                    mlp_thread_manager.wait_completion();
                    TIMING_ACCUMULATE(wait_for_mlp_finish_time, wait_for_mlp_t);

                    TIMING_START(part_single_node_t);
                    partition_single_node(partition_config, global_node_id, cur_line);
                    TIMING_ACCUMULATE(part_single_node_time, part_single_node_t);

                    buffer->update_neighbours_priority(cur_line);
                }


                // If buffer is full, partition either the top node or a batch of top nodes using MLP
                if (buffer->size() >= partition_config.max_pq_size) {
                    if (use_mlp) {

                        if (partition_config.parallel_mlp) {

                            // // Wait for MLP to finish
                            // TIMING_START(wait_for_mlp_t);
                            // mlp_thread_manager.wait_completion();
                            // TIMING_ACCUMULATE(wait_for_mlp_finish_time, wait_for_mlp_t);

                            // // CURRENTLY NOT WORKING
                            // buffer->loadTopNodesToBatch(batch_nodes, partition_config.stream_buffer_len, cur_batch_id);

                            // mlp_thread_manager.execute(partition_config, batch_nodes, cur_batch_id);

                            // cur_batch_id = partition_config.batch_manager->acquire_id();
                            // cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);
                        } else {
                            if (current_batch == nullptr) {
                                current_batch = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>();
                                current_batch->reserve(MIN(partition_config.stream_buffer_len, buffer->size()));
                            }

                            bool perform_mlp = false;
                            if (partition_config.batch_extraction_strategy == BATCH_EXTRACTION_STRATEGY_ALWAYS_TOP_NODE) {

                                LongNodeID top_node_id = buffer->deleteMax();
                                std::vector<LongNodeID> top_node_adj = std::move(buffer->get_adjacents(top_node_id));

                                if (partition_config.sep_batch_marker) {
                                    (*partition_config.stream_nodes_batch_marker)[top_node_id - 1] = cur_batch_marker;
                                } else {
                                    (*partition_config.stream_nodes_assign)[top_node_id - 1] = cur_batch_marker;
                                }

                                buffer->update_neighbours_priority(top_node_adj, false);
                                buffer->completely_remove_node(top_node_id);

                                current_batch->emplace_back(top_node_id, std::move(top_node_adj));
                                if (current_batch->size() >= partition_config.stream_buffer_len) {
                                    perform_mlp = true;
                                }
                            } else {
                                // Extract complete batch of nodes right away
                                buffer->loadTopNodesToBatch(current_batch, partition_config.stream_buffer_len, cur_batch_id);
                                perform_mlp = true;
                            }

                            if (perform_mlp) {
                                partition_config.nmbNodes = current_batch->size();
                                TIMING_START(t_mlp);
                                perform_mlp_on_batch(partition_config, current_batch, cur_batch_id);
                                TIMING_ACCUMULATE(mlp_time, t_mlp);

                                partition_config.batch_manager->release_id(cur_batch_id);
                                cur_batch_id = partition_config.batch_manager->acquire_id();
                                cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);
                            }

                        }

                    } else {
                        buffer->partitionTopNode();
                    }
                }


            } else {
                // debug_print("Restreaming: Adding node " + std::to_string(global_node_id) + " with degree " + std::to_string(degree));
                if (current_batch == nullptr) {
                    current_batch = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>();
                    current_batch->reserve(MIN(partition_config.stream_buffer_len, partition_config.remaining_stream_nodes+1));
                }
                current_batch->emplace_back(global_node_id, cur_line);
                if (current_batch->size() >= partition_config.stream_buffer_len || partition_config.remaining_stream_nodes == 0) {
                    partition_config.nmbNodes = current_batch->size();
                    // debug_print("[perform_mlp_on_batch] Restreaming: Processing batch of size " + std::to_string(current_batch->size()) + " with id " + std::to_string(cur_batch_id));

                    perform_mlp_on_batch(partition_config, current_batch, cur_batch_id);
                    // current_batch = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>();

                    partition_config.batch_manager->release_id(cur_batch_id);
                    cur_batch_id = partition_config.batch_manager->acquire_id();
                    cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);
                }
            }

        }
        cur_line.clear();
        TIMING_ACCUMULATE(first_phase_time, first_phase_t);
        (*partition_config.stream_in).close();

        if (partition_config.parallel_mlp) {
            // Wait for MLP to finish
            TIMING_START(wait_for_mlp_t);
            mlp_thread_manager.wait_completion();
            TIMING_ACCUMULATE(wait_for_mlp_finish_time, wait_for_mlp_t);

            mlp_time = mlp_thread_manager.get_mlp_time();
        }



        TIMING_START(second_phase_t);
        if (partition_config.restream_number == 0) {
            if ( use_mlp ) {
                while (!buffer->isEmpty()) {
                    if (current_batch == nullptr) {
                        current_batch = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>();
                        current_batch->reserve(MIN(partition_config.stream_buffer_len, buffer->size()));
                    }
                    buffer->loadTopNodesToBatch(current_batch, partition_config.stream_buffer_len, cur_batch_id);
                    TIMING_START(t_mlp);
                    partition_config.nmbNodes = current_batch->size();
                    perform_mlp_on_batch(partition_config, current_batch, cur_batch_id);
                    TIMING_ACCUMULATE(mlp_time, t_mlp);
                    // cur_batch_id = partition_config.batch_manager->acquire_id();
                    // cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);
                }
            } else {
                while (!buffer->isEmpty()) {
                    buffer->partitionTopNode();
                }
            }

            // Delete buffer
            updating_adj_time = buffer->get_update_adj_time();
            delete buffer;
        }

        if (partition_config.parallel_mlp) {
            // Wait for MLP to finish
            TIMING_START(wait_for_mlp_t);
            mlp_thread_manager.wait_completion();
            TIMING_ACCUMULATE(wait_for_mlp_finish_time, wait_for_mlp_t);

            mlp_time = mlp_thread_manager.get_mlp_time();
        }
        TIMING_ACCUMULATE(second_phase_time, second_phase_t);

        // if (partition_config.ghostkey_to_edges != NULL) {
        //     delete partition_config.ghostkey_to_edges;
        // }
        // if (partition_config.ghostkey_to_node != NULL) {
        //     delete partition_config.ghostkey_to_node;
        // }

        // debug_print("Finished restreaming pass: " + std::to_string(partition_config.restream_number));

    } // End of restreaming loop

    delete partition_config.batch_manager;


    double total_time = processing_t.elapsed();
    long maxRSS = getMaxRSS();


    if (partition_config.print_times) {
#ifdef ENABLE_TIME_MEASUREMENTS
        double sum_detailed;

        std::cout << "┌─────────────────────────┬───────────────┬───────────────┐" << std::endl;
        std::cout << "│ Metric                  │ Time (s)      │ Percentage    │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ Total time              │ " << std::setw(13) << std::fixed << std::setprecision(3) << total_time << " │ " << std::setw(13) << "100%" << " │" << std::endl;
        std::cout << "│ First phase time        │ " << std::setw(13) << std::fixed << std::setprecision(3) << first_phase_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (first_phase_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ Second phase time       │ " << std::setw(13) << std::fixed << std::setprecision(3) << second_phase_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (second_phase_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ IO time                 │ " << std::setw(13) << std::fixed << std::setprecision(3) << io_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (io_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ Buffer add node time    │ " << std::setw(13) << std::fixed << std::setprecision(3) << buffer_add_node_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (buffer_add_node_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ Updating adj time       │ " << std::setw(13) << std::fixed << std::setprecision(3) << updating_adj_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (updating_adj_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ Part single node time   │ " << std::setw(13) << std::fixed << std::setprecision(3) << part_single_node_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (part_single_node_time / total_time * 100) << "%" << " │" << std::endl;
        if (partition_config.parallel_mlp) {
            std::cout << "│ Wait for MLP finish     │ " << std::setw(13) << std::fixed << std::setprecision(3) << wait_for_mlp_finish_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (wait_for_mlp_finish_time / total_time * 100) << "%" << " │" << std::endl;
            std::cout << "│ Thread2: MLP time       │ " << std::setw(13) << std::fixed << std::setprecision(3) << mlp_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (mlp_time / total_time * 100) << "%" << " │" << std::endl;
            sum_detailed = io_time + buffer_add_node_time + updating_adj_time + wait_for_mlp_finish_time + part_single_node_time;
        } else {
            std::cout << "│ MLP time                │ " << std::setw(13) << std::fixed << std::setprecision(3) << mlp_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (mlp_time / total_time * 100) << "%" << " │" << std::endl;
            sum_detailed = io_time + buffer_add_node_time + updating_adj_time + mlp_time + part_single_node_time;
        }
        std::cout << "│ Sum of detailed times   │ " << std::setw(13) << std::fixed << std::setprecision(3) << sum_detailed << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (sum_detailed / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "└─────────────────────────┴───────────────┴───────────────┘" << std::endl;
#else
        std::cout << "Timing disabled - compile with -DENABLE_TIME_MEASUREMENTS to see detailed timing information" << std::endl;
#endif
    }
    FlatBufferWriter fb_writer;

    // Check if all nodes are assigned
    LongNodeID count_not_assigned = 0;
    for (LongNodeID i = 0; i < partition_config.number_of_nodes; i++) {
        if ((*partition_config.stream_nodes_assign)[i] > partition_config.k - 1) {
            count_not_assigned++;
            std::cout << "Node " << i+1 << " was not assigned to any partition. : " <<(*partition_config.stream_nodes_assign)[i] << std::endl;
        }
        // ASSERT_TRUE((*partition_config.stream_nodes_assign)[i] < TO_BE_PARTITIONED-10000);
    }
    if (count_not_assigned > 0) {
        std::cout << "Total nodes not assigned : " << count_not_assigned << std::endl;
    }

    graph_io_stream::streamEvaluatePartition(partition_config, graph_filename, total_edge_cut);
    double total_imbalance = qm.balance_full_stream(*partition_config.stream_blocks_weight);

    fb_writer.updateVertexPartitionResults(total_edge_cut, total_imbalance);


    double total_time_rounded = std::round(total_time * 1000.0) / 1000.0;
    std::cout << std::fixed << std::setprecision(3) << total_time_rounded
                << " " << maxRSS
                << " " << total_edge_cut
                << " " << std::defaultfloat << total_edge_cut / (double) partition_config.total_edges
                << " " << total_imbalance << std::endl;

    // write the partition to the disc
    std::stringstream filename;
    if (!partition_config.filename_output.compare("")) {
        filename << "tmppartition" << partition_config.k;
    } else {
        filename << partition_config.filename_output;
    }

    if (!partition_config.suppress_output) {
        graph_io_stream::writePartitionStream(partition_config, filename.str());
    } else {
        // std::cout << "No partition will be written as output." << std::endl;
    }

    if (partition_config.ghostkey_to_edges != NULL) {
        delete partition_config.ghostkey_to_edges;
    }

    if (partition_config.add_blocks_weight != NULL) {
        delete partition_config.add_blocks_weight;
    }
    if (partition_config.stream_nodes_assign != NULL) {
        delete partition_config.stream_nodes_assign;
    }
    if (partition_config.stream_nodes_batch_marker != NULL) {
        delete partition_config.stream_nodes_batch_marker;
    }
    if (partition_config.stream_blocks_weight != NULL) {
        delete partition_config.stream_blocks_weight;
    }
    if (partition_config.stream_in != NULL) {
        delete partition_config.stream_in;
    }

    fb_writer.updateResourceConsumption(buffer_io_time, model_construction_time, global_mapping_time, global_mapping_time, total_time, maxRSS);
    fb_writer.write(graph_filename, partition_config);

    return 0;
}




long getMaxRSS() {
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // The maximum resident set size is in kilobytes
        return usage.ru_maxrss;
    } else {
        std::cerr << "Error getting resource usage information." << std::endl;
        // Return a sentinel value or handle the error in an appropriate way
        return -1;
    }
}

// Function to extract the base filename without path and extension
std::string extractBaseFilename(const std::string &fullPath) {
    size_t lastSlash = fullPath.find_last_of('/');
    size_t lastDot = fullPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
        // Found a slash, extract the substring after the last slash
        return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    } else {
        // No slash found, just extract the substring before the last dot
        return fullPath.substr(0, lastDot);
    }
}



