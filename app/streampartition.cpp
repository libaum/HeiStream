// current best: ex2_v1.2
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


long getMaxRSS();

std::string extractBaseFilename(const std::string &fullPath);

void write_node_part_order_to_file(std::vector<std::string> &node_part_order) {
    std::ofstream order_file("hs_node_part_order");
    if (!order_file.is_open()) {
        std::cerr << "Failed to open node_part_order file for writing" << std::endl;
        return;
    }

    for (const auto& entry : node_part_order) {
        order_file << entry << "\n";
    }
    order_file.close();
}

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
    double global_mapping_time = 0;
    double buffer_io_time = 0;
    double io_time = 0;
    double buffer_add_node_time = 0;
    double model_construction_time = 0;
    double first_phase_time = 0;
    double second_phase_time = 0;
    double updating_adj_time = 0;
    double part_single_node_time = 0;
    double mlp_time = 0;
    double wait_for_mlp_finish_time = 0;

    PartitionConfig partition_config;
    std::string graph_filename;
    EdgeWeight total_edge_cut = 0;

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



    if (partition_config.write_node_part_order) {
        partition_config.node_part_order = new std::vector<std::string>();
    } else {
        partition_config.node_part_order = nullptr;
    }

    int &passes = partition_config.num_streams_passes;
    partition_config.count_misc1 = 0;
    partition_config.count_misc2 = 0;
    for (partition_config.restream_number = 0;
        partition_config.restream_number < passes; partition_config.restream_number++) {

        MLPThreadManager mlp_thread_manager;

        // ***************************** IO operations ***************************************
        io_t.restart();
        graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);

        double avg_block_size = static_cast<double>(partition_config.number_of_nodes) / partition_config.k;
        partition_config.max_block_weight = static_cast<int>(std::ceil((1.0 + partition_config.imbalance / 100) * avg_block_size));
        io_time += io_t.elapsed();

        Buffer buffer(partition_config, partition_config.max_pq_size);

        std::unique_ptr<buffered_input> ss2 = nullptr;
        std::vector<LongNodeID> cur_line;

        auto lines = std::make_unique<std::vector<std::string>>(1);

        // bool useFirstPhaseBuffer = partition_config.first_phase_buffer_len != 1;
        first_phase_t.restart();
        bool use_mlp = partition_config.stream_buffer_len != 1;
        while (partition_config.remaining_stream_nodes) {

            io_t.restart();
            // Load a line from the stream
            std::getline(*(partition_config.stream_in), (*lines)[0]);
            if ((*lines)[0][0] == '%') { // a comment in the file
                continue;
            }
            partition_config.remaining_stream_nodes--;
            LongNodeID global_node_id = ++partition_config.total_nodes_loaded;

            assert(global_node_id <= partition_config.number_of_nodes);

            ss2 = std::make_unique<buffered_input>(lines.get());
            ss2->simple_scan_line(cur_line);
            io_time += io_t.elapsed();

            unsigned degree = cur_line.size();

            if (degree >= partition_config.d_max || degree == 0) { //
                // Partition node directly if degree is too high or 0

                if (degree == 0) {
                    // Put node into partition with lowest weight
                    PartitionID best_partition = 0;
                    LongNodeID min_weight = partition_config.max_block_weight;

                    // Wait for MLP to finish
                    wait_for_mlp_t.restart();
                    mlp_thread_manager.wait_completion();
                    wait_for_mlp_finish_time += wait_for_mlp_t.elapsed();

                    for (PartitionID i = 0; i < partition_config.k; i++) {
                        if ((*partition_config.stream_blocks_weight)[i] < min_weight) {
                            min_weight = (*partition_config.stream_blocks_weight)[i];
                            best_partition = i;
                        }
                    }
                    (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;
                    (*partition_config.stream_blocks_weight)[best_partition]++;


                } else {

                    if (partition_config.write_node_part_order) {
                        std::string entry = std::to_string(global_node_id) + " " + std::to_string(-1) + " -> ";
                        (*partition_config.node_part_order).push_back(entry);
                    }

                    // Partition node directly
                    partition_config.count_misc1++;

                    // Wait for MLP to finish
                    wait_for_mlp_t.restart();
                    mlp_thread_manager.wait_completion();
                    wait_for_mlp_finish_time += wait_for_mlp_t.elapsed();

                    part_single_node_t.restart();
                    partition_single_node(partition_config, global_node_id, cur_line, true);
                    part_single_node_time += part_single_node_t.elapsed();

                    // Update neighbors
                    buffer.update_neighbours_priority(cur_line);
                }
                continue;

            } else if (buffer.size() >= partition_config.max_pq_size) {
                // Make space by removing node from queue by popping
                if (use_mlp) {

                    if (partition_config.parallel_mlp) {

                        // Wait for MLP to finish
                        wait_for_mlp_t.restart();
                        mlp_thread_manager.wait_completion();
                        wait_for_mlp_finish_time += wait_for_mlp_t.elapsed();

                        buffer.loadTopNodesToBatch(batch_nodes, partition_config.stream_buffer_len);

                        mlp_thread_manager.execute(partition_config, batch_nodes);
                    } else {
                        buffer.loadTopNodesToBatch(batch_nodes, partition_config.stream_buffer_len);

                        t_mlp.restart();
                        perform_mlp_on_batch(partition_config, batch_nodes);
                        mlp_time += t_mlp.elapsed();
                    }

                } else {
                    buffer.partitionTopNode();
                }
            }

            buffer_add_node_t.restart();
            buffer.addNode(global_node_id, cur_line);
            buffer_add_node_time += buffer_add_node_t.elapsed();

        }
        cur_line.clear();
        first_phase_time += first_phase_t.elapsed();
        (*partition_config.stream_in).close();


        second_phase_t.restart();
        if ( use_mlp ) {
            while (!buffer.isEmpty()) {

                if (partition_config.parallel_mlp) {

                    // Wait for MLP to finish
                    wait_for_mlp_t.restart();
                    mlp_thread_manager.wait_completion();
                    wait_for_mlp_finish_time += wait_for_mlp_t.elapsed();

                    buffer.loadTopNodesToBatch(batch_nodes, partition_config.stream_buffer_len);
                    mlp_thread_manager.execute(partition_config, batch_nodes);
                } else {
                    buffer.loadTopNodesToBatch(batch_nodes, partition_config.stream_buffer_len);
                    t_mlp.restart();
                    perform_mlp_on_batch(partition_config, batch_nodes);
                    mlp_time += t_mlp.elapsed();
                }


            }
        } else {
            while (!buffer.isEmpty()) {
                buffer.partitionTopNode();
            }
        }
        if (partition_config.parallel_mlp) {
            // Wait for MLP to finish
            wait_for_mlp_t.restart();
            mlp_thread_manager.wait_completion();
            wait_for_mlp_finish_time += wait_for_mlp_t.elapsed();

            mlp_time = mlp_thread_manager.get_mlp_time();
        }
        second_phase_time += second_phase_t.elapsed();
        updating_adj_time = buffer.get_update_adj_time();

    }
    double total_time = processing_t.elapsed();
    long maxRSS = getMaxRSS();


    if (partition_config.print_times) {
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
    }
    FlatBufferWriter fb_writer;

    // Check if all nodes are assigned
    for (LongNodeID i = 0; i < partition_config.number_of_nodes; i++) {
        ASSERT_TRUE((*partition_config.stream_nodes_assign)[i] < TO_BE_PARTITIONED);
    }

    graph_io_stream::streamEvaluatePartition(partition_config, graph_filename, total_edge_cut);
    fb_writer.updateVertexPartitionResults(total_edge_cut, qm.balance_full_stream(*partition_config.stream_blocks_weight));


    if (partition_config.write_node_part_order) {
        write_node_part_order_to_file(*partition_config.node_part_order);
    }

    double total_time_rounded = std::round(total_time * 100.0) / 100.0;
    std::cout << total_time_rounded << " " << total_edge_cut << " " << maxRSS << std::endl;

    // write the partition to the disc
    std::stringstream filename;
    if (!partition_config.filename_output.compare("")) {
        filename << "tmppartition" << partition_config.k;
    } else {
        filename << partition_config.filename_output;
    }

    if (!partition_config.suppress_output) {
        // graph_io_stream::writePartitionStream(partition_config, filename.str());
    } else {
        std::cout << "No partition will be written as output." << std::endl;
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



