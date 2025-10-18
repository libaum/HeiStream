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

#include "readerwriterqueue.h"
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

    timer processing_t, io_t, t_mlp, first_phase_t, second_phase_t, part_single_node_t, buffer_add_node_t;
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

    TIMING_DECLARE(double io_thread_time) = 0;
    TIMING_DECLARE(double pq_thread_time) = 0;
    TIMING_DECLARE(double partition_thread_time) = 0;

    // Wait / blocking times
    TIMING_DECLARE(double io_wait_time) = 0.0;
    TIMING_DECLARE(double pq_wait_time) = 0.0;
    TIMING_DECLARE(double partition_wait_time) = 0.0;

    // Timing variables for threads
    TIMING_DECLARE(double thread_io_time) = 0.0;
    TIMING_DECLARE(double thread_buffer_add_node_time) = 0.0;
    TIMING_DECLARE(double thread_part_single_node_time) = 0.0;
    TIMING_DECLARE(double thread_mlp_time) = 0.0;

    // Throughput metrics
    std::atomic<size_t> nodes_processed_io{0};
    std::atomic<size_t> nodes_processed_pq{0};
    std::atomic<size_t> tasks_processed_partition{0};

    // Monitor queue sizes
    std::atomic<size_t> max_input_queue_size{0};
    std::atomic<size_t> max_partition_queue_size{0};

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

    int &passes = partition_config.num_streams_passes;

    partition_config.batch_manager = new BatchIDManager(partition_config.max_active_batches);

    // Thread-safe queues and synchronization
    moodycamel::ReaderWriterQueue<ParsedLine>* input_queue_fast;
    moodycamel::ReaderWriterQueue<PartitionTask>* partition_queue_fast;
    if (!partition_config.alt_thread_queue) {
        input_queue_fast = new moodycamel::ReaderWriterQueue<ParsedLine>(partition_config.max_input_q_size);
        partition_queue_fast = new moodycamel::ReaderWriterQueue<PartitionTask>(100);
    }
    std::queue<ParsedLine> input_queue;
    std::queue<PartitionTask> partition_queue;
    std::mutex input_mutex, partition_mutex;
    std::condition_variable input_cv, partition_cv;

    std::atomic<bool> io_finished{false};
    std::atomic<bool> pq_finished{false};

    if  (partition_config.bb_ratio != UNDEFINED_BB_RATIO) {
        partition_config.batch_size = partition_config.max_buffer_size / partition_config.bb_ratio;
    }

    Buffer* buffer = new Buffer(partition_config, partition_config.max_buffer_size);

    auto push_to_partition_queue = [&](PartitionTask &&task) {

        if (!partition_config.alt_thread_queue) {
            bool succeeded = false;
            while (succeeded == false) {
                succeeded = partition_queue_fast->try_enqueue(std::move(task));  // Will only succeed if the queue has an empty slot (never allocates)
            }
        } else {
            {
                std::lock_guard<std::mutex> lock(partition_mutex);
                partition_queue.push(std::move(task));
            }
            partition_cv.notify_one();

        }


    };

    if (buffer != nullptr) {
        buffer->set_push_task_callback(push_to_partition_queue);
    }


    // Helper function to create single node partition task
    auto create_single_node_task = [&](LongNodeID node_id, std::vector<LongNodeID>&& adjacents) {
        if (partition_config.restream_number == 0) {
            // Check if node_id is currently active as ghost neighbor
            bool is_active_ghost_neighbor = false;
            if (partition_config.ghost_neighbors_enabled) {
                is_active_ghost_neighbor = partition_config.k <= (*partition_config.stream_nodes_assign)[node_id - 1]
                                                && (*partition_config.stream_nodes_assign)[node_id - 1] < 2 * partition_config.k;
            }
            buffer->update_neighbours_priority_parallel(
                adjacents,
                true,
                is_active_ghost_neighbor
            );
            (*partition_config.stream_nodes_assign)[node_id - 1] = TO_BE_PARTITIONED;
        }

        PartitionTask task(-1, {{node_id, std::move(adjacents)}});

        push_to_partition_queue(std::move(task));

    };


    for (partition_config.restream_number = 0; partition_config.restream_number < passes; partition_config.restream_number++) {


        // ***************************** IO operations ***************************************
        TIMING_START(io_t);
        graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);

        double avg_block_size = static_cast<double>(partition_config.number_of_nodes) / partition_config.k;
        partition_config.max_block_weight = static_cast<int>(std::ceil((1.0 + partition_config.imbalance / 100) * avg_block_size));
        TIMING_ACCUMULATE(io_time, io_t);

        bool use_mlp = partition_config.batch_size != 1;
        io_finished = false;
        pq_finished = false;

        // THREAD 1: IOReader
        std::thread io_reader([&]() {
            timer thread_timer, local_io_t;
            thread_timer.restart();

            std::vector<std::string> lines(1);
            buffered_input ss2(&lines);
            std::vector<LongNodeID> adjacents;
            adjacents.reserve(1000);

            partition_config.total_nodes_loaded = 0;
            ParsedLine parsed_line;

            timer wait_timer;
            while (partition_config.remaining_stream_nodes) {
                TIMING_START(local_io_t);

                // Load a line from the stream
                std::getline(*(partition_config.stream_in), lines[0]);
                if (lines[0][0] == '%') { // skip comments in the file
                    continue;
                }

                partition_config.remaining_stream_nodes--;
                LongNodeID global_node_id = ++partition_config.total_nodes_loaded;

                assert(global_node_id <= partition_config.number_of_nodes);

                ss2.simple_scan_line_fast(adjacents);
                TIMING_ACCUMULATE(thread_io_time, local_io_t);

                if ( partition_config.restream_number && !partition_config.restream_include_high_degree_nodes ) {
                    if ( adjacents.size() == 0 ||  adjacents.size() >= partition_config.d_max ) {
                        // In restreaming if not including high degree nodes, we do not process nodes with degree 0 or degree > d_max
                        adjacents = std::vector<LongNodeID>();
                    }
                }

                parsed_line = ParsedLine{global_node_id, adjacents};

                TIMING_START(wait_timer);
                if (!partition_config.alt_thread_queue) {
                    // Try enqueue in input_queue_fast
                    bool succeeded = false;
                    while (succeeded == false) {
                        succeeded = input_queue_fast->try_enqueue(std::move(parsed_line));  // Will only succeed if the queue has an empty slot (never allocates)
                    }

                } else {
                    {
                        std::unique_lock<std::mutex> lock(input_mutex);

                        // Wait if queue is too full (backpressure)
                        input_cv.wait(lock, [&] {
                            return input_queue.size() < partition_config.max_input_q_size;
                        });

                        TIMING_ACCUMULATE(io_wait_time, wait_timer);

                        // Monitor queue size
                        TIMING_MAX_UPDATE(max_input_queue_size, input_queue.size());

                        input_queue.push(std::move(parsed_line));
                    }
                    input_cv.notify_one();
                }


                TIMING_INCREMENT(nodes_processed_io);
                adjacents.clear();
            }

            // Update thread runtime
            TIMING_ACCUMULATE(io_thread_time, thread_timer);

            if (partition_config.alt_thread_queue) {
                io_finished = true;
                input_cv.notify_all();
            }
            (*partition_config.stream_in).close();

        });

        // THREAD 2: PQHandler
        std::thread pq_handler([&]() {
            timer thread_timer, local_buffer_t;

            // Batch management variables
            std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> current_batch;
            current_batch.reserve(partition_config.batch_size);

            size_t cur_batch_id = partition_config.batch_manager->acquire_id();
            PartitionID cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);

            LongNodeID total_nodes_processed = 0;

            // Helper function for batch creation from manually collected nodes
            auto create_batch_task = [&]() {
                // Create batch task with manually collected nodes
                PartitionTask task(cur_batch_id, std::move(current_batch));

                current_batch.clear(); // Reset for next batch
                current_batch.reserve(partition_config.batch_size);

                push_to_partition_queue(std::move(task));

                // Update batch id and marker for next batch
                cur_batch_id = partition_config.batch_manager->acquire_id();
                cur_batch_marker = partition_config.batch_manager->get_batch_marker(cur_batch_id);
            };



            // Helper function to manually extract top node from buffer
            auto extract_top_node_from_buffer = [&]() -> std::pair<LongNodeID, std::vector<LongNodeID>> {
                LongNodeID node_id = buffer->deleteMax();
                std::vector<LongNodeID> adjacents = std::move(buffer->get_adjacents(node_id));

                // Check if node_id is currently active as ghost neighbor
                bool is_active_ghost_neighbor = false;
                if (partition_config.ghost_neighbors_enabled) {
                    is_active_ghost_neighbor = partition_config.k <= (*partition_config.stream_nodes_assign)[node_id - 1]
                                                    && (*partition_config.stream_nodes_assign)[node_id - 1] < 2 * partition_config.k;
                }

                buffer->update_neighbours_priority_parallel(adjacents, false, is_active_ghost_neighbor);
                if (partition_config.sep_batch_marker) {
                    (*partition_config.stream_nodes_batch_marker)[node_id - 1] = cur_batch_marker;
                } else {
                    (*partition_config.stream_nodes_assign)[node_id - 1] = cur_batch_marker;
                }

                buffer->completely_remove_node(node_id);
                return {node_id, std::move(adjacents)};
            };



            // Helper function to try extracting one node for current batch
            auto add_node_to_batch = [&](LongNodeID& node_id, std::vector<LongNodeID>&& node_adjacents) -> bool {
                // auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                current_batch.emplace_back(node_id, std::move(node_adjacents));

                // Check if batch is complete
                if (current_batch.size() >= partition_config.batch_size) {
                    create_batch_task();
                    return true; // Batch completed
                }
                return false; // Batch not yet complete or no extraction possible
            };

            timer wait_timer;
            while (true) {
                if (total_nodes_processed >= partition_config.number_of_nodes) {
                    break;
                }
                // Process incoming nodes from input queue
                ParsedLine parsed_line;

                if (!partition_config.alt_thread_queue) {
                    int removed = input_queue_fast->try_dequeue(parsed_line);
                    if (!removed) continue;
                } else {
                    // Wait for input or termination
                    TIMING_START(wait_timer);
                    {
                        std::unique_lock<std::mutex> lock(input_mutex);

                        input_cv.wait(lock, [&] { return !input_queue.empty() || io_finished; });

                        TIMING_ACCUMULATE(pq_wait_time, wait_timer);

                        if (input_queue.empty() && io_finished) {
                            break;
                        }

                        parsed_line = std::move(input_queue.front());
                        input_queue.pop();
                        input_cv.notify_one();

                    }
                }
                total_nodes_processed++;

                TIMING_START(thread_timer);

                // Process the input node if we have one
                LongNodeID global_node_id = parsed_line.node_id;
                std::vector<LongNodeID>& adjacents = parsed_line.neighbors;
                unsigned degree = adjacents.size();

                if (partition_config.restream_number) {
                    if (degree == 0) { // In case of restreaming this means degree 0 or degree > d_max
                        create_single_node_task(global_node_id, std::vector<LongNodeID>());
                        TIMING_INCREMENT(nodes_processed_pq);
                        continue; // Skip empty nodes
                    }

                    // In restreaming, we do not add nodes to the buffer, but directly create partition tasks
                    add_node_to_batch(global_node_id, std::move(adjacents));

                } else {
                    if (degree >= partition_config.d_max || degree == 0) {
                        create_single_node_task(global_node_id, std::move(adjacents));
                        TIMING_ACCUMULATE(pq_thread_time, thread_timer);
                        TIMING_INCREMENT(nodes_processed_pq);
                        continue;
                    }

                    TIMING_START(local_buffer_t);
                    bool added_to_buffer = buffer->addNode(global_node_id, adjacents);
                    TIMING_ACCUMULATE(thread_buffer_add_node_time, local_buffer_t);

                    if (!added_to_buffer) {
                        create_single_node_task(global_node_id, std::move(adjacents));

                    }

                    if (buffer->size() > partition_config.max_buffer_size) {
                        if (use_mlp) {
                            if (partition_config.batch_extraction_strategy != BATCH_EXTRACTION_STRATEGY_ALWAYS_TOP_NODE) {
                                // Extract top nodes and their neighbors from the buffer and create batch task
                                buffer->loadTopNodesAndNeighborsToBatch(current_batch, partition_config.batch_size, cur_batch_id);
                                create_batch_task();

                            } else {
                                auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                                add_node_to_batch(node_id, std::move(node_adjacents));
                            }

                        } else {
                            // Handle non-MLP mode buffer overflow
                            auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                            create_single_node_task(node_id, std::move(node_adjacents));
                        }

                    }
                }
                TIMING_ACCUMULATE(pq_thread_time, thread_timer);
                TIMING_INCREMENT(nodes_processed_pq);
            }

            TIMING_START(thread_timer);


            // Handle remaining buffer
            if (partition_config.restream_number == 0) {
                if (use_mlp) {
                    // Process remaining nodes in batches using interleaved approach
                    while (!buffer->isEmpty()) {
                        if (partition_config.batch_extraction_strategy != BATCH_EXTRACTION_STRATEGY_ALWAYS_TOP_NODE) {
                            // Extract top nodes and their neighbors from the buffer and create batch task
                            buffer->loadTopNodesAndNeighborsToBatch(current_batch, partition_config.batch_size, cur_batch_id);
                            create_batch_task();

                        } else {
                            auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                            add_node_to_batch(node_id, std::move(node_adjacents));
                        }
                    }
                } else {
                    // Process remaining nodes individually
                    while (!buffer->isEmpty()) {
                        auto [node_id, node_adjacents] = extract_top_node_from_buffer();
                        create_single_node_task(node_id, std::move(node_adjacents));
                    }
                }
            }

            // Complete any remaining batch
            if (!current_batch.empty()) {
                create_batch_task();
            }

            TIMING_ACCUMULATE(pq_thread_time, thread_timer);

            if (partition_config.alt_thread_queue) {
                pq_finished = true;
                partition_cv.notify_all();
            }
        });

        // THREAD 3: PartitionWorker
        std::thread partition_worker([&]() {
            timer thread_timer, local_part_t, local_mlp_t;

            LongNodeID total_nodes_partitioned = 0;

            timer wait_timer;
            while (true) {
                if (total_nodes_partitioned >= partition_config.number_of_nodes) {
                    break; // All nodes processed
                }

                PartitionTask task;

                if (!partition_config.alt_thread_queue) {
                    int removed = partition_queue_fast->try_dequeue(task);
                    if (!removed) continue;

                } else {
                    TIMING_START(wait_timer);

                    {
                        std::unique_lock<std::mutex> lock(partition_mutex);
                        partition_cv.wait(lock, [&] { return !partition_queue.empty() || pq_finished; });

                        TIMING_ACCUMULATE(partition_wait_time, wait_timer);

                        if (partition_queue.empty() && pq_finished) {
                            break;
                        }

                        TIMING_MAX_UPDATE(max_partition_queue_size, partition_queue.size());

                        task = std::move(partition_queue.front());
                        partition_queue.pop();
                    }

                }


                // TIMING_START(thread_timer);

                if (task.batch_id == -1) {
                    // Single node partitioning
                    LongNodeID node_id = task.nodes[0].first;
                    auto& adjacents = task.nodes[0].second;
                    total_nodes_partitioned++;

                    if (partition_config.restream_number) {
                        TIMING_ACCUMULATE(thread_part_single_node_time, local_part_t);
                        continue; // Skip processing single nodes in restreaming mode
                    }

                    // Single node partitioning
                    TIMING_START(local_part_t);

                    if (adjacents.size() == 0) {
                        PartitionID best_partition = 0;
                        LongNodeID min_weight = partition_config.max_block_weight;

                        for (PartitionID i = 0; i < partition_config.k; i++) {
                            if ((*partition_config.stream_blocks_weight)[i] < min_weight) {
                                min_weight = (*partition_config.stream_blocks_weight)[i];
                                best_partition = i;
                            }
                        }
                        (*partition_config.stream_nodes_assign)[node_id - 1] = best_partition;
                        (*partition_config.stream_blocks_weight)[best_partition]++;
                    } else {
                        partition_single_node(partition_config, node_id, adjacents);
                    }

                    TIMING_ACCUMULATE(thread_part_single_node_time, local_part_t);

                } else {
                    // Batch partitioning (batch_id >= 0)
                    TIMING_START(local_mlp_t);

                    // Convert our manual batch format to the expected format
                    auto batch_ptr = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>(std::move(task.nodes));
                    partition_config.nmbNodes = batch_ptr->size();
                    total_nodes_partitioned += batch_ptr->size();
                    if (partition_config.restream_number) {
                        // log_memory_usage("Before perform_mlp_on_batch for batch " + std::to_string(task.batch_id));
                    }
                    perform_mlp_on_batch(partition_config, batch_ptr, task.batch_id);
                    partition_config.batch_manager->release_id(task.batch_id);

                    TIMING_ACCUMULATE(thread_mlp_time, local_mlp_t);

                    task.nodes.clear();
                }

                TIMING_ACCUMULATE(partition_thread_time, thread_timer);
                TIMING_INCREMENT(tasks_processed_partition);
            }
        });

        // Join all threads
        io_reader.join();
        pq_handler.join();
        partition_worker.join();


        if (partition_config.restream_number == 0) {
            // graph_io_stream::writePartitionStream(partition_config, "tmp1"); //, filename.str());
            if (buffer != nullptr) {
                updating_adj_time = buffer->get_update_adj_time();
                delete buffer;
                buffer = nullptr;
            }
        }

        if (partition_config.ghostkey_to_edges != NULL) {
            delete partition_config.ghostkey_to_edges;
            partition_config.ghostkey_to_edges = NULL;
        }

        if (partition_config.ghostkey_to_node != NULL) {
            delete partition_config.ghostkey_to_node;
            partition_config.ghostkey_to_node = NULL;
        }

        if (partition_config.ghostglobal_to_ghostkey != NULL) {
            delete partition_config.ghostglobal_to_ghostkey;
            partition_config.ghostglobal_to_ghostkey = NULL;
        }

        if (partition_config.stream_in != NULL) {
            delete partition_config.stream_in;
            partition_config.stream_in = NULL;
        }

    }

    if (partition_config.batch_manager != nullptr) {
        delete partition_config.batch_manager;
    }

    // Update timing variables (thread-safe, since all threads have finished)
    TIMING_ADD(io_time, thread_io_time);
    TIMING_ADD(buffer_add_node_time, thread_buffer_add_node_time);
    TIMING_ADD(part_single_node_time, thread_part_single_node_time);
    TIMING_ADD(mlp_time, thread_mlp_time);

    TIMING_ACCUMULATE(first_phase_time, first_phase_t);
    TIMING_ACCUMULATE(second_phase_time, second_phase_t);


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
        std::cout << "│ MLP time                │ " << std::setw(13) << std::fixed << std::setprecision(3) << mlp_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (mlp_time / total_time * 100) << "%" << " │" << std::endl;
        sum_detailed = io_time + buffer_add_node_time + updating_adj_time + mlp_time + part_single_node_time;
        std::cout << "│ Sum of detailed times   │ " << std::setw(13) << std::fixed << std::setprecision(3) << sum_detailed << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (sum_detailed / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ THREAD RUNTIME TRACKING │               │               │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ IOReader thread time    │ " << std::setw(13) << std::fixed << std::setprecision(3) << io_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (io_thread_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ PQHandler thread time   │ " << std::setw(13) << std::fixed << std::setprecision(3) << pq_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (pq_thread_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ PartitionWorker time    │ " << std::setw(13) << std::fixed << std::setprecision(3) << partition_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (partition_thread_time / total_time * 100) << "%" << " │" << std::endl;
        double total_thread_time = io_thread_time + pq_thread_time + partition_thread_time;
        std::cout << "│ Total thread time       │ " << std::setw(13) << std::fixed << std::setprecision(3) << total_thread_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (total_thread_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ BOTTLENECK ANALYSIS     │               │               │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ IO wait time            │ " << std::setw(13) << std::fixed << std::setprecision(3) << io_wait_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (io_wait_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ PQ wait time            │ " << std::setw(13) << std::fixed << std::setprecision(3) << pq_wait_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (pq_wait_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "│ Partition wait time     │ " << std::setw(13) << std::fixed << std::setprecision(3) << partition_wait_time << " │ " << std::setw(12) << std::fixed << std::setprecision(0) << (partition_wait_time / total_time * 100) << "%" << " │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ Nodes/sec IO            │ " << std::setw(13) << std::fixed << std::setprecision(0) << (TIMING_LOAD(nodes_processed_io) / total_time) << " │               │" << std::endl;
        std::cout << "│ Nodes/sec PQ            │ " << std::setw(13) << std::fixed << std::setprecision(0) << (TIMING_LOAD(nodes_processed_pq) / total_time) << " │               │" << std::endl;
        std::cout << "│ Tasks/sec Partition     │ " << std::setw(13) << std::fixed << std::setprecision(0) << (TIMING_LOAD(tasks_processed_partition) / total_time) << " │               │" << std::endl;
        std::cout << "├─────────────────────────┼───────────────┼───────────────┤" << std::endl;
        std::cout << "│ Max input queue size    │ " << std::setw(13) << std::fixed << std::setprecision(0) << TIMING_LOAD(max_input_queue_size) << " │               │" << std::endl;
        std::cout << "│ Max partition queue size│ " << std::setw(13) << std::fixed << std::setprecision(0) << TIMING_LOAD(max_partition_queue_size) << " │               │" << std::endl;
        std::cout << "└─────────────────────────┴───────────────┴───────────────┘" << std::endl;
#else
        std::cout << "Timing disabled - compile with -DENABLE_TIME_MEASUREMENTS to see detailed timing information" << std::endl;
#endif
    }
    FlatBufferWriter fb_writer;

    // Check if all nodes are assigned
    for (LongNodeID i = 0; i < partition_config.number_of_nodes; i++) {
        ASSERT_TRUE((*partition_config.stream_nodes_assign)[i] < TO_BE_PARTITIONED-10000);
    }

    graph_io_stream::streamEvaluatePartition(partition_config, graph_filename, total_edge_cut);
    double total_imbalance = qm.balance_full_stream(*partition_config.stream_blocks_weight);
    fb_writer.updateVertexPartitionResults(total_edge_cut, total_imbalance);

    double total_time_rounded = std::round(total_time * 1000.0) / 1000.0;
    std::cout << std::fixed << std::setprecision(3) << total_time_rounded
                << " " << maxRSS
                << " " << total_edge_cut
                << " " << std::defaultfloat << total_edge_cut / (double) partition_config.total_edges
                << std::endl;

    // write the partition to the disc
    if (!partition_config.suppress_output) {
        graph_io_stream::writePartitionStream(partition_config);
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



