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
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <atomic>


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

std::string extractBaseFilename(const std::string &fullPath);

void write_new_npo_entry(PartitionConfig& partition_config, LongNodeID global_node_id) {
    if (partition_config.write_node_part_order) {
        std::string entry = std::to_string(global_node_id) + " " + std::to_string(partition_config.k) + " -> ";
        (*partition_config.node_part_order).push_back(entry);
    }
}

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

    double io_thread_time = 0;
    double pq_thread_time = 0;
    double partition_thread_time = 0;

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
    }

    int &passes = partition_config.num_streams_passes;
    partition_config.count_misc1 = 0;
    partition_config.count_misc2 = 0;

    for (partition_config.restream_number = 0;
        partition_config.restream_number < passes; partition_config.restream_number++) {

        // MLPThreadManager mlp_thread_manager;
        partition_config.batch_manager = new BatchIDManager(partition_config.max_active_batches);

        // ***************************** IO operations ***************************************
        io_t.restart();
        graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);

        double avg_block_size = static_cast<double>(partition_config.number_of_nodes) / partition_config.k;
        partition_config.max_block_weight = static_cast<int>(std::ceil((1.0 + partition_config.imbalance / 100) * avg_block_size));
        io_time += io_t.elapsed();

        Buffer buffer(partition_config, partition_config.max_pq_size);

        // Thread-safe queues and synchronization
        std::queue<ParsedLine> input_queue;
        std::queue<PartitionTask> partition_queue;
        std::mutex input_mutex, partition_mutex;
        std::condition_variable input_cv, partition_cv;
        std::atomic<bool> io_finished{false};
        std::atomic<bool> pq_finished{false};

        // Timing variables für threads
        double thread_io_time = 0.0;
        double thread_buffer_add_node_time = 0.0;
        double thread_part_single_node_time = 0.0;
        double thread_mlp_time = 0.0;


        // std::atomic<double> thread_wait_for_mlp_time{0.0};

        bool use_mlp = partition_config.stream_buffer_len != 1;

        // THREAD 1: IOReader
        std::thread io_reader([&]() {
            timer thread_timer, local_io_t;
            thread_timer.restart();

            std::vector<std::string> lines(1);
            buffered_input ss2(&lines);
            std::vector<LongNodeID> adjacents;
            adjacents.reserve(1000);

            while (partition_config.remaining_stream_nodes) {
                local_io_t.restart();

                // Load a line from the stream
                std::getline(*(partition_config.stream_in), lines[0]);
                if (lines[0][0] == '%') { // skip comments in the file
                    continue;
                }

                partition_config.remaining_stream_nodes--;
                LongNodeID global_node_id = ++partition_config.total_nodes_loaded;

                assert(global_node_id <= partition_config.number_of_nodes);

                ss2.simple_scan_line_fast(adjacents);

                thread_io_time += local_io_t.elapsed();

                // Create ParsedLine and push to queue with backpressure
                ParsedLine parsed_line{global_node_id, adjacents};
                {
                    std::unique_lock<std::mutex> lock(input_mutex);

                    // Wait if queue is too full (backpressure)
                    input_cv.wait(lock, [&] {
                        return input_queue.size() < partition_config.max_input_q_size;
                    });

                    input_queue.push(std::move(parsed_line));
                }
                input_cv.notify_one();

                adjacents.clear();
            }

            // Update thread runtime
            io_thread_time += thread_timer.elapsed();

            io_finished = true;
            input_cv.notify_all();
            (*partition_config.stream_in).close();
        });

        // THREAD 2: PQHandler
        std::thread pq_handler([&]() {
            timer thread_timer, local_buffer_t;

            // Helper function for batch creation
            auto create_batch_task = [&]() {
                size_t batch_id = partition_config.batch_manager->acquire_id(); // Necessary as a unique marker that a node is contained in a specific batch in partition_config->stream_nodes_assign

                std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> *batch_nodes_ptr;
                buffer.loadTopNodesToBatch(batch_nodes_ptr, partition_config.stream_buffer_len, batch_id);

                PartitionTask task(batch_id, batch_nodes_ptr);

                {
                    std::lock_guard<std::mutex> lock(partition_mutex);
                    partition_queue.push(std::move(task));
                }
                partition_cv.notify_one();
            };

            // Helper function to create single node partition task
            auto create_single_node_task = [&](LongNodeID node_id, std::vector<LongNodeID>&& adjacents) {
                // TODO: handle case with part_adj_directly correctly, currently should not work
                // - create queue for adjacents for neighbors that should be partitioned directly,
                //    -> then after finishing this function, parition them in while loop?
                // NO: actually just add a task for them -> i will handle this in the PartitionWorker

                // buffer.update_neighbours_priority(adjacents, true);
                buffer.update_neighbours_priority(
                    adjacents,
                    true,
                    &partition_queue,
                    &partition_mutex,
                    &partition_cv,
                    partition_config.stream_nodes_assign
                );
                (*partition_config.stream_nodes_assign)[node_id - 1] = TO_BE_PARTITIONED;

                PartitionTask task(
                    -1,
                    std::vector<BatchNode>{{node_id, std::move(adjacents)}}
                );

                {
                    std::lock_guard<std::mutex> lock(partition_mutex);
                    partition_queue.push(std::move(task));
                }
                partition_cv.notify_one();
            };

            // Helper function for task creation for the top node in the buffer
            auto create_task_for_top_node = [&]() {
                LongNodeID node_id = buffer.deleteMax();
                // Get the adjacents of the node to be partitioned, move to a new vector
                std::vector<LongNodeID> adjacents_copy = std::move(buffer.get_adjacents(node_id));
                buffer.completely_remove_node(node_id);

                create_single_node_task(node_id, std::move(adjacents_copy));
            };

            while (true) {
                ParsedLine parsed_line;

                // Wait for input or termination
                {
                    std::unique_lock<std::mutex> lock(input_mutex);
                    input_cv.wait(lock, [&] { return !input_queue.empty() || io_finished; });

                    if (input_queue.empty() && io_finished) {
                        break;
                    }

                    parsed_line = std::move(input_queue.front());
                    input_queue.pop();

                    // Notify IOReader that queue has space (backpressure relief)
                    input_cv.notify_one();
                }
                thread_timer.restart();

                LongNodeID global_node_id = parsed_line.node_id;
                std::vector<LongNodeID>& adjacents = parsed_line.neighbors;
                unsigned degree = adjacents.size();

                if (degree >= partition_config.d_max || degree == 0) {
                    create_single_node_task(global_node_id, std::move(adjacents));
                    pq_thread_time += thread_timer.elapsed();
                    continue;
                }

                // Check if new node has a higher buffer score than max score in buffer
                local_buffer_t.restart();
                bool added_to_buffer = buffer.addNode(global_node_id, adjacents);

                thread_buffer_add_node_time += local_buffer_t.elapsed();

                if (!added_to_buffer) {
                    // Means that the buffer is full and the current node has the highest score -> create partition task for it
                    create_single_node_task(global_node_id, std::move(adjacents));
                }

                // If buffer is full, create batch task
                if (buffer.size() > partition_config.max_pq_size) {
                    if (use_mlp) {
                        create_batch_task();
                    } else {
                        create_task_for_top_node();
                    }
                }
                pq_thread_time += thread_timer.elapsed();
            }

            thread_timer.restart();
            // Handle remaining buffer
            if (use_mlp) {
                while (!buffer.isEmpty()) {
                    create_batch_task();
                }
            } else {
                while (!buffer.isEmpty()) {
                    create_task_for_top_node();
                }
            }
            // Update thread runtime
            pq_thread_time += thread_timer.elapsed();

            pq_finished = true;
            partition_cv.notify_all();
        });

        // THREAD 3: PartitionWorker
        std::thread partition_worker([&]() {
            timer thread_timer, local_part_t, local_mlp_t;

            while (true) {
                PartitionTask task;

                // Wait for partition task or termination
                {
                    std::unique_lock<std::mutex> lock(partition_mutex);
                    partition_cv.wait(lock, [&] { return !partition_queue.empty() || pq_finished; });

                    if (partition_queue.empty() && pq_finished) {
                        break;
                    }

                    task = std::move(partition_queue.front());
                    partition_queue.pop();
                }

                thread_timer.restart();

                if (task.batch_id == -1) {
                    // Single node partitioning (stack-allocated)
                    local_part_t.restart();
                    auto& node_data = task.nodes[0]; // node_data is a pair<LongNodeID, std::vector<LongNodeID>>

                    if (node_data.second.size() == 0) {
                        // Put node into partition with lowest weight
                        PartitionID best_partition = 0;
                        LongNodeID min_weight = partition_config.max_block_weight;

                        for (PartitionID i = 0; i < partition_config.k; i++) {
                            if ((*partition_config.stream_blocks_weight)[i] < min_weight) {
                                min_weight = (*partition_config.stream_blocks_weight)[i];
                                best_partition = i;
                            }
                        }
                        (*partition_config.stream_nodes_assign)[node_data.first - 1] = best_partition;
                        (*partition_config.stream_blocks_weight)[best_partition]++;
                    } else {
                        partition_single_node(partition_config, node_data.first, node_data.second);
                    }

                    thread_part_single_node_time += local_part_t.elapsed();

                } else {
                    // Batch partitioning (heap-allocated)
                    local_mlp_t.restart();

                    auto batch_ptr = task.heap_batch_nodes;
                    partition_config.nmbNodes = batch_ptr->size();
                    perform_mlp_on_batch(partition_config, batch_ptr, task.batch_id);
                    partition_config.batch_manager->release_id(task.batch_id);

                    thread_mlp_time += local_mlp_t.elapsed();
                }
                // Update thread runtime
                partition_thread_time += thread_timer.elapsed();
            }

        });

        // Join all threads
        io_reader.join();
        pq_handler.join();
        partition_worker.join();

        // Update timing variables (thread-sicher, da alle Threads beendet sind)
        io_time += thread_io_time;
        buffer_add_node_time += thread_buffer_add_node_time;
        part_single_node_time += thread_part_single_node_time;
        mlp_time += thread_mlp_time;

        first_phase_time += first_phase_t.elapsed();
        second_phase_time += second_phase_t.elapsed();
        updating_adj_time = buffer.get_update_adj_time();

        if (partition_config.print_times) {
            buffer.print_pq_statistics();
        }

        delete partition_config.batch_manager;
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
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ THREAD RUNTIME TRACKING │               │               │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ IOReader thread time    │ " << std::setw(13) << std::fixed << std::setprecision(3) << io_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (io_thread_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ PQHandler thread time   │ " << std::setw(13) << std::fixed << std::setprecision(3) << pq_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (pq_thread_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ PartitionWorker time    │ " << std::setw(13) << std::fixed << std::setprecision(3) << partition_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (partition_thread_time / total_time * 100) << "%" << " │" << std::endl;
        double total_thread_time = io_thread_time + pq_thread_time + partition_thread_time;
        std::cout << "│ Total thread time       │ " << std::setw(13) << std::fixed << std::setprecision(3) << total_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (total_thread_time / total_time * 100) << "%" << " │" << std::endl;
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

    double total_time_rounded = std::round(total_time * 1000.0) / 1000.0;
    std::cout << std::fixed << std::setprecision(3) << total_time_rounded;
    std::cout << " " << maxRSS;
    std::cout << " " << total_edge_cut;
    std::cout << " " << std::defaultfloat << total_edge_cut / (double) partition_config.total_edges << std::endl;

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



