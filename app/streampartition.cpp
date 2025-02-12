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

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
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

#define MIN(A, B) ((A)<(B))?(A):(B)
#define MAX(A, B) ((A)>(B))?(A):(B)

int cnt_part_adj_directly = 0;

// CUTTANA HYPERPARAMETERS
const double THETA = 2;
const int D_MAX = 1000;
const bool PARTITION_ADJ_DIRECTLY_ENABLED = true;
const bool HEISTREAM_INTEGRATION_ENABLED = true;

// Constants  for bucket queue
const float MAX_BUFFER_SCORE = 3.0f;
const LongNodeID INVALID_NODE = std::numeric_limits<LongNodeID>::max();

class buffer_item {
public:
    float buffer_score;
    std::vector<LongNodeID> line;
    int num_adj_partitioned;

    buffer_item(std::vector<LongNodeID> line, float buffer_score, int num_adj_partitioned)
        : line(line), buffer_score(buffer_score), num_adj_partitioned(num_adj_partitioned) {}

    void clear() {
        line.clear();
        line.shrink_to_fit();
        buffer_score = UNDEFINED_LONGNODE;
        num_adj_partitioned = 0;
    }
};

static int discretize_score(PartitionConfig &partition_config, float score) {
    // Use round instead of floor to handle precision better
    return static_cast<int>(std::round(score * partition_config.bq_disc_factor));
}


void loadTopNodesToBatch(PartitionConfig &partition_config, bucket_pq &pq, std::vector<buffer_item> &node_id_to_buffer_item, std::vector<std::vector<LongNodeID>> *&input, LongNodeID batch_size);
void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<std::vector<LongNodeID>> *&input);

void config_multibfs_initial_partitioning(PartitionConfig &partition_config);
void partition_node(PartitionConfig &partition_config, std::vector<LongNodeID> &line);

float calc_buffer_score(PartitionConfig &partition_config,
                        std::vector<LongNodeID> &cur_line,
                        int &cnt_adj_partitioned) {
    int global_node_id = cur_line[0];
    float degree = cur_line.size() - 1;
    if (degree == 0) {
        return 0;
    }
    // float cnt_adj_partitioned = 0;
    ASSERT_TRUE(cnt_adj_partitioned == 0);
    bool adj_is_partitioned;
    for (LongNodeID &global_adj_id : cur_line) {
        if (global_adj_id == global_node_id)
            continue;
        adj_is_partitioned = (*partition_config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION;
        if (adj_is_partitioned) {
            cnt_adj_partitioned++;
        }
    }
    float buffer_score = degree / D_MAX + THETA * (float)cnt_adj_partitioned / degree; // Range: [0, 3] (first term: [0, 1], second term: [0, THETA])
    return buffer_score;
}

void assert_neighbors_partitioned(PartitionConfig &partition_config, std::vector<LongNodeID> &line, bool all_should_be_partitioned) {
    bool one_not_partitioned = false;
    for (LongNodeID &global_adj_id : line) {
        if (global_adj_id == line[0])
            continue;
        if ((*partition_config.stream_nodes_assign)[global_adj_id - 1] == INVALID_PARTITION) {
            if (all_should_be_partitioned) {
                ASSERT_TRUE(false);
            }

            one_not_partitioned = true;
            break;
        }
    }
    if (!all_should_be_partitioned) {
        ASSERT_TRUE(one_not_partitioned);
    }
}


// Update the priority value of the neighbours of the node that was just partitioned in the priority queue
void update_neighbours_priority(PartitionConfig &partition_config,
                                std::vector<LongNodeID> &line,
                                std::vector<buffer_item> &node_id_to_buffer_item,
                                bucket_pq &pq, bool part_adj_directly=PARTITION_ADJ_DIRECTLY_ENABLED) {

    // Ensure that the line contains at least the node id (and possibly neighbors).
    if( line.size() == 0 || line.size() == 1) {
        // std::cout << "Line is empty." << std::endl;
        return;
    }

    for (auto it = line.begin() + 1; it != line.end(); ++it) {
        LongNodeID adj_id = *it;
        // Check if the neighbour is already partitioned (if MLP batch partitioning, node is only virutally partitioned and therefore line is empty)
        // bool is_partitioned = ((*partition_config.stream_nodes_assign)[adj_id - 1] != INVALID_PARTITION)
        //     || ((*partition_config.node_in_current_block)[adj_id - 1] == 1);

        // if (pq.contains(adj_id)) {
        bool is_partitioned = (*partition_config.stream_nodes_assign)[adj_id - 1] != INVALID_PARTITION; // || node_id_to_buffer_item[adj_id - 1].line.size() == 0;
        auto &adj_buffer_item = node_id_to_buffer_item[adj_id - 1];
        if (!is_partitioned && adj_buffer_item.buffer_score != UNDEFINED_LONGNODE) { // pq.contains(adj_id)) { // {
            if (adj_buffer_item.line.size() < 2) {
                continue;
            }

            adj_buffer_item.num_adj_partitioned++;
            int adj_degree = adj_buffer_item.line.size() - 1;

            // Check if all neighbours of the neighbour are partitioned, if so, partition the neighbour
            if (part_adj_directly && adj_degree > 3 && adj_degree == adj_buffer_item.num_adj_partitioned) {
                // assert_neighbors_partitioned(partition_config, adj_buffer_item.line, true);
                cnt_part_adj_directly++;

                // std::cout << "Partitioning node " << adj_id << " directly, buffer score: " << adj_buffer_item.buffer_score + THETA / adj_degree << ", degree: " << adj_degree << std::endl;
                // for (LongNodeID &adj_global_id : adj_buffer_item.line) {
                //     std::cout << adj_global_id << " ";            }
                // std::cout << std::endl;
                pq.deleteNode(adj_id);
                partition_node(partition_config, adj_buffer_item.line);

                // Update neighbors and clear buffer item
                update_neighbours_priority(partition_config, adj_buffer_item.line, node_id_to_buffer_item, pq);
                adj_buffer_item.clear();
            } else {
                // Update buffer score of neighbours
                // assert_neighbors_partitioned(partition_config, adj_buffer_item.line, false);
                float updated_buffer_score = adj_buffer_item.buffer_score + THETA / adj_degree;
                pq.increaseKey(adj_id, discretize_score(partition_config, updated_buffer_score));
                adj_buffer_item.buffer_score = updated_buffer_score;
            }
        }
    }
}

void partition_node(PartitionConfig &partition_config, std::vector<LongNodeID> &line) {
    LongNodeID global_node_id = line[0];
    std::vector<int> hash_map(partition_config.k, 0);

    for (auto it = line.begin() + 1; it != line.end(); ++it) {
        LongNodeID adj_id = *it;
        PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id - 1];
        if (adj_part != INVALID_PARTITION) {
            hash_map[adj_part]++;
        }
    }

    // Step 3: Iterate over partitions to compute FENNEL scores
    PartitionID best_partition = 0;
    float best_score = std::numeric_limits<float>::lowest();
    bool feasible_partition_found = false;

    for (PartitionID cur_partition = 0; cur_partition < partition_config.k; ++cur_partition) {
        NodeWeight current_block_weight = (*partition_config.stream_blocks_weight)[cur_partition];
        // Skip or penalize partitions that are already "full"
        if (current_block_weight >= partition_config.max_block_weight) {
            // This partition is not feasible anymore
            continue;
            // Or set score = -infinity, but "continue" is simpler
        }

        float edge_gain = hash_map[cur_partition];
        NodeWeight partition_load = (*partition_config.stream_blocks_weight)[cur_partition];
        float load_penalty = partition_config.fennel_alpha_gamma * random_functions::approx_sqrt(partition_load);
        float score = edge_gain - load_penalty;

        if (score > best_score) {
            best_score = score;
            best_partition = cur_partition;
            feasible_partition_found = true;
        }
    }

    if (!feasible_partition_found) {
        // If no feasible partition is found, assign to the partition with the least load
        std::cout << "No feasible partition found for node " << global_node_id << ". Assigning to partition with least load." << std::endl;
        best_partition = std::min_element(partition_config.stream_blocks_weight->begin(), partition_config.stream_blocks_weight->end()) - partition_config.stream_blocks_weight->begin();
    } else {
        // Step 4: Assign the node to the partition with the best score
        (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;

        // Step 5: Update partition load
        (*partition_config.stream_blocks_weight)[best_partition]++;
    }
}

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
    std::vector <std::vector<LongNodeID>> *input = NULL;

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

    timer first_pass_t, second_pass_t, updating_adj_t, partitioning_t, calc_buffer_score_t;
    double first_pass_time = 0;
    double second_pass_time = 0;
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
        bucket_pq pq(static_cast<int>(std::floor(MAX_BUFFER_SCORE * partition_config.bq_disc_factor)) + 1);

        std::vector<buffer_item> node_id_to_buffer_item(partition_config.number_of_nodes, buffer_item(std::vector<LongNodeID>(0), UNDEFINED_LONGEDGE, 0));

        LongNodeID node_counter = 0;
        std::unique_ptr<buffered_input> ss2 = nullptr;
        std::vector<LongNodeID> cur_line;

        auto lines = std::make_unique<std::vector<std::string>>(1);

        first_pass_t.restart();
        while (partition_config.remaining_stream_nodes) {
            // Load a line from the stream
            std::getline(*(partition_config.stream_in), (*lines)[0]);
            if ((*lines)[0][0] == '%') { // a comment in the file
                continue;
            }
            partition_config.remaining_stream_nodes--;
            LongNodeID global_node_id = ++partition_config.total_nodes_loaded;

            ASSERT_TRUE(global_node_id <= partition_config.number_of_nodes);

            ss2 = std::make_unique<buffered_input>(lines.get());
            ss2->simple_scan_line(cur_line);
            cur_line.insert(cur_line.begin(), global_node_id);

            int degree = cur_line.size() - 1;
            if (degree > D_MAX || degree == 0) {
                // Partition node directly
                partitioning_t.restart();
                partition_node(partition_config, cur_line);
                partitioning_time += partitioning_t.elapsed();
                updating_adj_t.restart();

                // Update neighbors and clear buffer item
                update_neighbours_priority(partition_config, cur_line, node_id_to_buffer_item, pq);
                updating_adj_time += updating_adj_t.elapsed();
                node_id_to_buffer_item[global_node_id - 1].clear();
                continue;
            } else if (pq.size() >= partition_config.max_pq_size) {
                // Make space by removing node from queue by popping
                if (partition_config.first_phase_buffer_len == 1) {
                    LongNodeID node_id_to_remove = pq.deleteMax();
                    ASSERT_TRUE((*partition_config.stream_nodes_assign)[node_id_to_remove - 1] == INVALID_PARTITION);

                    // Partition the node
                    partitioning_t.restart();
                    partition_node(partition_config, node_id_to_buffer_item[node_id_to_remove - 1].line);
                    partitioning_time += partitioning_t.elapsed();
                    updating_adj_t.restart();

                    // Update neighbors and clear buffer item
                    update_neighbours_priority(partition_config, node_id_to_buffer_item[node_id_to_remove - 1].line, node_id_to_buffer_item, pq);
                    updating_adj_time += updating_adj_t.elapsed();
                    node_id_to_buffer_item[node_id_to_remove - 1].clear();
                } else {
                    loadTopNodesToBatch(partition_config, pq, node_id_to_buffer_item, input, partition_config.first_phase_buffer_len);
                    perform_mlp_on_batch(partition_config, input);
                }
            }

            // Calculate priority of the node and push into BucketQueue
            calc_buffer_score_t.restart();

            float buffer_score = calc_buffer_score(partition_config, cur_line, node_id_to_buffer_item[global_node_id - 1].num_adj_partitioned);

            calc_buffer_score_time += calc_buffer_score_t.elapsed();
            pq.insert(global_node_id, discretize_score(partition_config, buffer_score));
            node_id_to_buffer_item[global_node_id - 1].line = cur_line;
            node_id_to_buffer_item[global_node_id - 1].buffer_score = buffer_score;
        }
        cur_line.clear();
        first_pass_time += first_pass_t.elapsed();
        (*partition_config.stream_in).close();

        second_pass_t.restart();
        while (!pq.empty()) {
            if (HEISTREAM_INTEGRATION_ENABLED) {
                loadTopNodesToBatch(partition_config, pq, node_id_to_buffer_item, input, partition_config.second_phase_buffer_len);
                // loadTopNodesAboveThresholdToBatch(partition_config, pq, node_id_to_buffer_item, input);
                perform_mlp_on_batch(partition_config, input);
            } else {
                LongNodeID node_id_to_partition = pq.deleteMax();
                ASSERT_TRUE((*partition_config.stream_nodes_assign)[node_id_to_partition - 1] == INVALID_PARTITION);

                bool is_already_partitioned = (*partition_config.stream_nodes_assign)[node_id_to_partition - 1] != INVALID_PARTITION;
                if (is_already_partitioned) {
                    continue;
                }

                // Partition the node
                partitioning_t.restart();
                auto &buffer_item_to_be_partitioned = node_id_to_buffer_item[node_id_to_partition - 1];
                partition_node(partition_config, buffer_item_to_be_partitioned.line);
                partitioning_time += partitioning_t.elapsed();

                // Update neighbors and clear buffer item
                updating_adj_t.restart();
                update_neighbours_priority(partition_config, buffer_item_to_be_partitioned.line, node_id_to_buffer_item, pq);
                updating_adj_time += updating_adj_t.elapsed();
                buffer_item_to_be_partitioned.clear();
            }
        }
        second_pass_time += second_pass_t.elapsed();
    }
    // std::cout << "First pass time: " << first_pass_time << std::endl;
    // std::cout << "Second pass time: " << second_pass_time << std::endl;
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
    if (partition_config.node_in_current_block != NULL) {
        delete partition_config.node_in_current_block;
    }
    if (partition_config.stream_nodes_assign != NULL) {
        delete partition_config.stream_nodes_assign;
    }
    if (partition_config.local_to_global_map != NULL) {
        delete partition_config.local_to_global_map;
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

// Select top partition_config.second_phase_buffer_len number of nodes from bucket queues and load into the batch (input)
void loadTopNodesToBatch(PartitionConfig &partition_config,
                       bucket_pq &pq,
                       std::vector<buffer_item> &node_id_to_buffer_item,
                       std::vector<std::vector<LongNodeID>> *&input,
                       LongNodeID batch_size) {

    // Initialize the partition configuration
    partition_config.nmbNodes = MIN(batch_size, pq.size());
    // std::fill(partition_config.node_in_current_block->begin(), partition_config.node_in_current_block->end(), 0);

    input = new std::vector<std::vector<LongNodeID>>(partition_config.nmbNodes);

    // 1. Extract the top batch_size nodes from the queue
    int node_counter = 0;
    // std::cout << "##################### -------- LOADING NEW BATCH -------- #####################" << std::endl;
    while (node_counter < partition_config.nmbNodes && !pq.empty()) {
        LongNodeID node_id = pq.deleteMax();
        ASSERT_TRUE((*partition_config.stream_nodes_assign)[node_id - 1] == INVALID_PARTITION);

        (*partition_config.node_in_current_block)[node_id - 1] = 1;

        auto &buffer_item = node_id_to_buffer_item[node_id - 1];

        (*input)[node_counter] = buffer_item.line;

        // Update neighbors and clear buffer item
        update_neighbours_priority(partition_config, buffer_item.line, node_id_to_buffer_item, pq, false);

        buffer_item.clear();

        node_counter++;
    }
    // std::cout << std::endl;
}

// A function to do multi-level partitioning the nodes in the batch (input)
void perform_mlp_on_batch(PartitionConfig &partition_config, std::vector<std::vector<LongNodeID>> *&input) {
    // Initialize the partition configuration
    graph_access *G = new graph_access();
    quality_metrics qm;
    balance_configuration bc;

    // ***************************** build model ***************************************
    G->set_partition_count(partition_config.k);
    graph_io_stream::createModel(partition_config, *G, input);
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

    delete G;
}
