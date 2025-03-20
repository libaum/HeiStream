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
#include <fstream>
#include <sstream>
#include <memory>
#include <optional>

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



long getMaxRSS();

// void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<std::vector<LongNodeID>> *&input);
void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<LongNodeID> *&input_idxs, Buffer &buffer);

void config_multibfs_initial_partitioning(PartitionConfig &partition_config);

void assert_neighbors_partitioned(PartitionConfig &partition_config, std::vector<LongNodeID> &line, bool all_should_be_partitioned) {
    bool one_not_partitioned = false;
    for (LongNodeID &global_adj_id : line) {
        if (global_adj_id == line[0])
            continue;
        if ((*partition_config.stream_nodes_assign)[global_adj_id - 1] == INVALID_PARTITION) {
            if (all_should_be_partitioned) {
                assert(false);
            }

            one_not_partitioned = true;
            break;
        }
    }
    if (!all_should_be_partitioned) {
        assert(one_not_partitioned);
    }
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
    PartitionConfig partition_config;
    std::string graph_filename;
    timer t, processing_t, io_t, model_t;
    EdgeWeight total_edge_cut = 0;
    double global_mapping_time = 0;
    double buffer_mapping_time = 0;
    double buffer_io_time = 0;
    double model_construction_time = 0;
    quality_metrics qm;
    balance_configuration bc;
    // std::vector <std::vector<LongNodeID>> *input = NULL;
    std::vector<LongNodeID> *input_idxs = NULL;


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

    timer first_phase_t, second_phase_t, updating_adj_t, partitioning_t, calc_buffer_score_t;
    double first_phase_time = 0;
    double second_phase_time = 0;
    double updating_adj_time = 0;
    double partitioning_time = 0;
    double calc_buffer_score_time = 0;

    int &passes = partition_config.num_streams_passes;
    for (partition_config.restream_number = 0;
        partition_config.restream_number < passes; partition_config.restream_number++) {

        // ***************************** IO operations ***************************************
        io_t.restart();
        graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);

        double avg_block_size = static_cast<double>(partition_config.number_of_nodes) / partition_config.k;
        partition_config.max_block_weight = static_cast<int>(std::ceil((1.0 + partition_config.imbalance / 100) * avg_block_size));

        buffer_io_time += io_t.elapsed();

        // bucket_pq pq(static_cast<int>(std::floor(MAX_BUFFER_SCORE * partition_config.bq_disc_factor)) + 1, partition_config.number_of_nodes);
        // std::vector<PQItem> node_id_to_buffer_item(partition_config.number_of_nodes); // ex2_v1
        // std::unordered_map<LongNodeID, PQItem> node_id_to_buffer_item2; // = std::unordered_map<LongNodeID, PQItem>(0);

        Buffer buffer(partition_config);

        LongNodeID node_counter = 0;
        std::unique_ptr<buffered_input> ss2 = nullptr;
        std::vector<LongNodeID> cur_line;

        auto lines = std::make_unique<std::vector<std::string>>(1);

        first_phase_t.restart();

        // bool useFirstPhaseBuffer = partition_config.first_phase_buffer_len != 1;
        bool use_mlp = partition_config.stream_buffer_len != 1;
        while (partition_config.remaining_stream_nodes) {

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

            int degree = cur_line.size();
            if (degree > D_MAX || degree == 0) {
                // Partition node directly if degree is too high or 0
                partitioning_t.restart();
                partition_single_node(partition_config, global_node_id, cur_line);
                partitioning_time += partitioning_t.elapsed();
                updating_adj_t.restart();

                // Update neighbors
                buffer.update_neighbours_priority(cur_line);
                updating_adj_time += updating_adj_t.elapsed();
                continue;
            } else if (buffer.size() >= partition_config.max_pq_size) {
                // Make space by removing node from queue by popping
                if (use_mlp) {
                    buffer.loadTopNodesToBatch(input_idxs, partition_config.stream_buffer_len);
                    perform_mlp_on_batch(partition_config, input_idxs, buffer);
                } else {
                    buffer.partitionTopNode();
                }
            }

            buffer.addNode(global_node_id, cur_line);
        }
        cur_line.clear();
        first_phase_time += first_phase_t.elapsed();
        (*partition_config.stream_in).close();


        second_phase_t.restart();
        // bool useSecondPhaseBuffer = partition_config.second_phase_buffer_len != 1;
        if ( use_mlp ) {
            while (!buffer.isEmpty()) {
                buffer.loadTopNodesToBatch(input_idxs, partition_config.stream_buffer_len);
                perform_mlp_on_batch(partition_config, input_idxs, buffer);
            }
        } else {
            while (!buffer.isEmpty()) {
                buffer.partitionTopNode();
            }
        }
        second_phase_time += second_phase_t.elapsed();
    }
    // std::cout << "First phase time: " << first_phase_time << std::endl;
    // std::cout << "Second phase time: " << second_phase_time << std::endl;
    // std::cout << "Updating adj time: " << updating_adj_time << std::endl;
    // std::cout << "Partitioning time: " << partitioning_time << std::endl;
    // std::cout << "Calc buffer score time: " << calc_buffer_score_time << std::endl;
    double total_time = processing_t.elapsed();
    long maxRSS = getMaxRSS();
    FlatBufferWriter fb_writer;

    // Check if all nodes are assigned
    for (LongNodeID i = 0; i < partition_config.number_of_nodes; i++) {
        ASSERT_TRUE((*partition_config.stream_nodes_assign)[i] != INVALID_PARTITION);
        // if ((*partition_config.stream_nodes_assign)[i] == INVALID_PARTITION) {
        // 	std::cout << "Node " << i << " is not assigned." << std::endl;
        // }
    }

    graph_io_stream::streamEvaluatePartition(partition_config, graph_filename, total_edge_cut);
    fb_writer.updateVertexPartitionResults(total_edge_cut, qm.balance_full_stream(*partition_config.stream_blocks_weight));

    double total_time_rounded = std::round(total_time * 100.0) / 100.0;
    std::cout << total_time_rounded << " " << total_edge_cut << " " << maxRSS << std::endl; // << " " << cnt_part_adj_directly << std::endl;

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
    if (partition_config.node_in_current_block != nullptr) {
        delete partition_config.node_in_current_block;
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


void config_multibfs_initial_partitioning(PartitionConfig &partition_config) {
    if (partition_config.initial_part_multi_bfs && partition_config.curr_batch >= 2) {
        partition_config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
    }
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



// A function to do multi-level partitioning the nodes in the batch (input)
void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<LongNodeID> *&input_idxs, Buffer &buffer) {
    // Initialize the partition configuration
    graph_access *G = new graph_access();
    quality_metrics qm;
    balance_configuration bc;

    // ***************************** build model ***************************************
    G->set_partition_count(partition_config.k);
    partition_config.local_to_global_map = new std::vector<NodeID>(partition_config.nmbNodes, 0);
    graph_io_stream::createModel(partition_config, *G, input_idxs, buffer);
    graph_io_stream::countAssignedNodes(partition_config);
    graph_io_stream::prescribeBufferInbalance(partition_config);
    bool already_fully_partitioned = (partition_config.restream_vcycle && partition_config.restream_number);
    bc.configurate_balance(partition_config, *G, already_fully_partitioned || !partition_config.stream_initial_bisections);
    config_multibfs_initial_partitioning(partition_config);

    // ***************************** perform partitioning ***************************************
    graph_partitioner partitioner;
    partitioner.perform_partitioning(partition_config, *G);

    // ***************************** permanent assignment ***************************************
    graph_io_stream::generalizeStreamPartition(partition_config, *G);

    delete partition_config.local_to_global_map;
    delete G;
}
