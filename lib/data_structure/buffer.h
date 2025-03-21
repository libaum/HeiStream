#ifndef BUFFER_NKJSAF9
#define BUFFER_NKJSAF9

#include <iostream>
#include <memory>
#include <optional>
#include <vector>

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "macros_assertions.h"
#include "partition/partition_config.h"
#include "random_functions.h"
#include <sparsehash/dense_hash_map>

#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#define MAX(A, B) ((A) > (B)) ? (A) : (B)

// CUTTANA HYPERPARAMETERS
const double THETA = 2;
const int D_MAX = 1000;
const bool PARTITION_ADJ_DIRECTLY_ENABLED = true;

// Constants  for bucket queue
const float MAX_BUFFER_SCORE = 1 + THETA;

inline void partition_single_node(PartitionConfig &partition_config, LongNodeID global_node_id, std::vector<LongNodeID> &line) {

    std::vector<int> hash_map(partition_config.k, 0);
    for (LongNodeID adj_id : line) {
        PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id - 1];
        if (adj_part != INVALID_PARTITION) {
            hash_map[adj_part]++;
        }
    }

    // Iterate over partitions to compute FENNEL scores
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
        // Assign the node to the partition with the best score
        (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;

        // Update partition load
        (*partition_config.stream_blocks_weight)[best_partition]++;
    }
}

class Buffer {
private:
    bucket_pq pq;
    // std::vector<PQItem> node_id_to_buffer_item;
    // google::dense_hash_map<LongNodeID, PQItem> node_id_to_buffer_item2;
    std::unordered_map<LongNodeID, PQItem> node_id_to_buffer_item2;
    // std::unordered_map<LongNodeID, PQItem, UInt64Hash> node_id_to_buffer_item2;
    // int cnt_part_adj_directly;
    PartitionConfig &config;

public:
    Buffer(PartitionConfig &partition_config)
        : pq(static_cast<int>(std::floor(MAX_BUFFER_SCORE * partition_config.bq_disc_factor)) + 1, partition_config.number_of_nodes),
          // node_id_to_buffer_item(partition_config.number_of_nodes),
          config(partition_config) {
        // node_id_to_buffer_item2.set_deleted_key("");
        // node_id_to_buffer_item2.set_empty_key(-1);

        size_t expected_size = partition_config.max_pq_size;
        float max_load_factor = 0.7f; // Standard ist 1.0
        size_t bucket_count = std::ceil(expected_size / max_load_factor);

        node_id_to_buffer_item2.reserve(expected_size);
        node_id_to_buffer_item2.max_load_factor(max_load_factor);
        // node_id_to_buffer_item2.rehash(bucket_count);
    }

    // cnt_part_adj_directly(0)

    ~Buffer() {
        // Cleanup - Delete line pointers
        // for (auto& item : node_id_to_buffer_item) {
        //     if (item.is_valid() && item.line != nullptr) {
        //         delete item.line;
        //         item.line = nullptr;
        //     }
        // }
    }

    int discretize_score(float score) {
        // Use round instead of floor to handle precision better
        return static_cast<int>(std::round(score * config.bq_disc_factor));
    }

    // Adds a node to the buffer
    void addNode(LongNodeID global_node_id, std::vector<LongNodeID> &line) {

        // float buffer_score = calc_buffer_score(global_node_id, line, node_id_to_buffer_item[global_node_id - 1].num_adj_partitioned);
        int num_adj_partitioned = 0;
        float buffer_score = calc_buffer_score(global_node_id, line, num_adj_partitioned);

        node_id_to_buffer_item2[global_node_id] = PQItem(buffer_score, line, num_adj_partitioned);
        // node_id_to_buffer_item2.emplace(
        //     global_node_id,
        //     PQItem(buffer_score, line, num_adj_partitioned)
        // );
        pq.insert(global_node_id, discretize_score(buffer_score));

        // node_id_to_buffer_item[global_node_id - 1].line = new std::vector<LongNodeID>(line);
        // node_id_to_buffer_item[global_node_id - 1].buffer_score = buffer_score;
    }

    // Removes and partitions the node with the highest priority
    void partitionTopNode() {
        LongNodeID node_id = pq.deleteMax();
        assert((*config.stream_nodes_assign)[node_id - 1] == INVALID_PARTITION);

        // Partition the node
        // auto &line = *node_id_to_buffer_item[node_id - 1].line;
        auto &line = *node_id_to_buffer_item2[node_id].line;
        // pq.getKey(node_id);
        partition_single_node(config, node_id, line);

        // Update neighbors and clear buffer item
        update_neighbours_priority(line);

        // node_id_to_buffer_item[node_id - 1].make_invalid();
        // node_id_to_buffer_item[node_id - 1].clean();
        // node_id_to_buffer_item2[node_id].make_invalid();
        removeNode(node_id);
    }

    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesToBatch(std::vector<LongNodeID> *&input_idxs, LongNodeID batch_size) {
        // Initialize the partition configuration
        config.nmbNodes = MIN(batch_size, pq.size());
        input_idxs = new std::vector<LongNodeID>(config.nmbNodes);

        // Extract the top batch_size number of nodes from the queue
        int node_counter = 0;
        while (node_counter < config.nmbNodes && !pq.empty()) {
            LongNodeID node_id = pq.deleteMax();
            assert((*config.stream_nodes_assign)[node_id - 1] == INVALID_PARTITION);

            (*config.node_in_current_block)[node_id - 1] = 1;

            auto &buffer_item = node_id_to_buffer_item2[node_id];
            // pq.getKey(node_id);
            //  auto &buffer_item = node_id_to_buffer_item[node_id - 1];

            // Update neighbors and make buffer item invalid
            update_neighbours_priority(*buffer_item.line, false);

            (*input_idxs)[node_counter] = node_id;
            // node_id_to_buffer_item[node_id - 1].make_invalid();
            buffer_item.make_invalid();
            node_counter++;
        }
    }

    // Berechnet den Score für einen Knoten
    float calc_buffer_score(LongNodeID global_node_id, std::vector<LongNodeID> &cur_line, int &cnt_adj_partitioned) {
        // int global_node_id = cur_line[0];
        float degree = cur_line.size();
        if (degree == 0) {
            return 0;
        }
        // float cnt_adj_partitioned = 0;
        assert(cnt_adj_partitioned == 0);
        bool adj_is_partitioned;
        for (LongNodeID &global_adj_id : cur_line) {
            adj_is_partitioned = (*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION;
            if (adj_is_partitioned) {
                cnt_adj_partitioned++;
            }
        }
        float buffer_score = degree / D_MAX + THETA * (float)cnt_adj_partitioned / degree; // Range: [0, 3] (first term: [0, 1], second term: [0, THETA])
        return buffer_score;
    }

    // Update the priority value of the neighbours of the node that was just partitioned in the priority queue
    void update_neighbours_priority(std::vector<LongNodeID> &line, bool part_adj_directly = PARTITION_ADJ_DIRECTLY_ENABLED) {
        // Ensure that the line contains at least one neighbor.
        if (line.size() == 0) {
            return;
        }

        for (LongNodeID adj_id : line) {
            bool is_partitioned = (*config.stream_nodes_assign)[adj_id - 1] != INVALID_PARTITION;

            auto it = node_id_to_buffer_item2.find(adj_id);
            // pq.contains(adj_id);

            // if (!is_partitioned && node_id_to_buffer_item[adj_id - 1].is_valid()) { // pq.contains(adj_id)) { // {
            if (!is_partitioned && it != node_id_to_buffer_item2.end() && it->second.is_valid()) {
                // if (!is_partitioned && node_id_to_buffer_item2.count(adj_id) && node_id_to_buffer_item2[adj_id].is_valid()) {
                // auto &adj_buffer_item = node_id_to_buffer_item[adj_id - 1];
                auto &adj_buffer_item = it->second; // node_id_to_buffer_item2[adj_id];
                // pq.getKey(adj_id);
                int adj_degree = adj_buffer_item.line->size();
                if (adj_degree == 0) {
                    continue;
                }

                adj_buffer_item.num_adj_partitioned++;

                // Check if all neighbours of the neighbour are partitioned, if so, partition the neighbour
                if (part_adj_directly && adj_degree > 3 && adj_degree == adj_buffer_item.num_adj_partitioned) {
                    // assert_neighbors_partitioned(config, adj_buffer_item.line, true);
                    // cnt_part_adj_directly++;
                    pq.deleteNode(adj_id);
                    partition_single_node(config, adj_id, *adj_buffer_item.line);

                    // Update neighbors and clear buffer item
                    update_neighbours_priority(*adj_buffer_item.line);

                    // node_id_to_buffer_item[adj_id - 1].make_invalid();
                    // node_id_to_buffer_item[adj_id - 1].clean();
                    // adj_buffer_item.make_invalid();
                    // node_id_to_buffer_item2[adj_id].clean();
                    removeNode(adj_id);
                } else {
                    // Update buffer score of neighbours
                    // assert_neighbors_partitioned(config, adj_buffer_item.line, false);
                    float updated_buffer_score = adj_buffer_item.buffer_score + THETA / adj_degree;
                    float score = discretize_score(updated_buffer_score);
                    if (score < 0) {
                        std::cout << "Negative score: " << score << std::endl;
                        std::cout << "Updated buffer score: " << updated_buffer_score << std::endl;
                        std::cout << "Degree: " << adj_degree << std::endl;
                        std::cout << "node id: " << adj_id << std::endl;
                        std::cout << "neighbors:" << adj_buffer_item.line->size() << std::endl;
                        std::cout << "contains: " << node_id_to_buffer_item2.count(adj_id) << std::endl;
                        std::cout << "valid: " << node_id_to_buffer_item2[adj_id].is_valid() << std::endl;
                    }
                    pq.increaseKey(adj_id, score);
                    adj_buffer_item.buffer_score = updated_buffer_score;
                }
            }
        }
    }

    // Helper-Methoden
    bool isEmpty() const { return pq.empty(); }
    size_t size() const { return pq.size(); }

    PQItem &getBufferItem(LongNodeID node_id) {
        // return node_id_to_buffer_item[node_id - 1];
        // pq.getKey(node_id);
        return node_id_to_buffer_item2[node_id];
    }

    std::vector<LongNodeID> &getLine(LongNodeID node_id) {
        // return *node_id_to_buffer_item[node_id - 1].line;
        // pq.getKey(node_id);
        return *node_id_to_buffer_item2[node_id].line;
    }

    void removeNode(LongNodeID node_id) {
        // node_id_to_buffer_item[node_id - 1].clean();
        // pq.getKey(node_id);
        node_id_to_buffer_item2.erase(node_id);
    }
};

#endif // BUFFER_NKJSAF9