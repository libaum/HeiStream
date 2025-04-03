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
// #include <sparsehash/dense_hash_map>

#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#define MAX(A, B) ((A) > (B)) ? (A) : (B)

// CUTTANA HYPERPARAMETERS
const double THETA = 2;
const int D_MAX = 1000;
const bool PARTITION_ADJ_DIRECTLY_ENABLED = true;


inline void partition_single_node(PartitionConfig &partition_config, LongNodeID global_node_id, std::vector<LongNodeID> &adjacents, bool is_hub = false) {

    partition_config.count_misc2++;
    std::vector<int> hash_map(partition_config.k, 0);
    int cnt_future_neighbors = 0;
    for (LongNodeID adj_id : adjacents) {
        PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id - 1];
        if (adj_part != INVALID_PARTITION) {
            hash_map[adj_part]++;
        } else {
            cnt_future_neighbors++;
        }
    }

    // Iterate over partitions to compute FENNEL scores
    PartitionID best_partition = 0;
    float best_score = std::numeric_limits<float>::lowest();
    bool feasible_partition_found = false;

    float f = 1; // Balance-Penalty-Faktor for Hubs
    if (is_hub) {
        // f -= ((float) cnt_future_neighbors /  adjacents.size());
        // std::cout << "Node " << global_node_id << " is a hub. f = " << f << std::endl;
    }
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
        float score = edge_gain - f * load_penalty;

        // if (is_hub) {
        //     std::cout << "Partition " << cur_partition << " has edge gain " << edge_gain << " and load penalty " << load_penalty << "partitionLoad: " << partition_load << ". Score: " << score << std::endl;
        // }
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
    }
    // Assign the node to the partition with the best score
    (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;

    // Update partition load
    (*partition_config.stream_blocks_weight)[best_partition]++;

}

// Effiziente Implementierung von pow für Fließkomma-Exponenten
inline float fast_pow(float base, float exponent) {
    return std::exp(exponent * std::log(base));
}

// Optimierte Version mit Spezialisierung für häufige beta-Werte
inline float fast_pow_specialized(float base, float exponent) {
    // Optimierungen für häufige Exponenten
    if (exponent == 2.0f) {
        return base * base;
    }
    else if (exponent == 1.0f) {
        return base;
    }
    // else if (exponent == 0.5f) {
    //     return random_functions::approx_sqrt(base);
    // }
    // else if (exponent == 0.0f) {
    //     return 1.0f;
    // }
    // Allgemeiner Fall für andere Exponenten
    return std::exp(exponent * std::log(base));
}

class Buffer {
private:
    bucket_pq pq;
    // std::vector<PQItem> node_id_to_buffer_item;
    // google::dense_hash_map<LongNodeID, PQItem> node_id_to_buffer_item2;
    // std::unordered_map<LongNodeID, PQItem> node_id_to_buffer_item2;
    PartitionConfig &config;

    LongNodeID total_degree_sum;
    LongNodeID node_counter;
    float current_beta;
    // std::unordered_map<LongNodeID, PQItem, UInt64Hash> node_id_to_buffer_item2;
    // int cnt_part_adj_directly;

public:
    Buffer(PartitionConfig &partition_config)
        :   config(partition_config),
            pq(static_cast<unsigned>(std::floor(get_max_buffer_score(partition_config) * partition_config.bq_disc_factor)) + 1,
                partition_config.number_of_nodes, partition_config.max_pq_size) {
          // node_id_to_buffer_item(partition_config.number_of_nodes),
        // node_id_to_buffer_item2.set_deleted_key("");
        // node_id_to_buffer_item2.set_empty_key(-1);

        size_t expected_size = config.max_pq_size;
        float max_load_factor = 0.7f; // Standard ist 1.0
        size_t bucket_count = std::ceil(expected_size / max_load_factor);

        // node_id_to_buffer_item2.reserve(expected_size);
        // node_id_to_buffer_item2.max_load_factor(max_load_factor);

        current_beta = config.haa_beta;
        total_degree_sum = 0;
        node_counter = 0;
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

    static float get_max_buffer_score(const PartitionConfig& cfg) {
        switch (cfg.buffer_score_type) {
            case BUFFER_SCORE_ANC:
                return 1;
            case BUFFER_SCORE_ANC2:
                return 1;
            case BUFFER_SCORE_HAA:
                return std::max(cfg.haa_theta, 1.0f);
            case BUFFER_SCORE_CMS:
                return 1;
            case BUFFER_SCORE_NSS:
                return 1;
            case BUFFER_SCORE_GTS:
                return 1;
            case BUFFER_SCORE_CBS:
            case BUFFER_SCORE_CBS2:
            default:
                return 3;
        }
    }

    // Berechnet den Score für einen Knoten
    float calc_buffer_score(LongNodeID global_node_id, const std::vector<LongNodeID> &adjacents, unsigned &cnt_adj_partitioned) {
        // int global_node_id = adjacents[0];
        if (adjacents.size() == 0) {
            return 0;
        }

        float degree = static_cast<float>(adjacents.size());
        float buffer_score;
        assert(cnt_adj_partitioned == 0);

        switch (config.buffer_score_type) {
            case BUFFER_SCORE_ANC: // Assigned Neighbor Count
            {
                for (const LongNodeID &global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    }
                }
                buffer_score = cnt_adj_partitioned /  degree;
                break;
            }
            case BUFFER_SCORE_ANC2: // Assigned Neighbor Count
            {
                float alpha = 0.6f;
                float beta = 1.5f;

                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    }
                }

                float neighbor_term = cnt_adj_partitioned / degree;  // avoid div by 0
                float deg_term = pow(degree / config.d_max, beta);

                buffer_score = alpha * neighbor_term + (1.0f - alpha) * deg_term;
                break;
            }
            case BUFFER_SCORE_HAA:
            {
                // Zähle, wie viele Nachbarn schon partitioniert sind.
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    }
                }
                // Berechne den adaptiven Gewichtungsfaktor r(v)
                float r = std::pow(degree / static_cast<float>(config.d_max), current_beta);
                // float r = fast_pow(degree / config.d_max, current_beta);
                // float r = fast_pow_specialized(degree / config.d_max, current_beta);

                // Berechne den Score: für niedrige Degree (r~0) dominiert der Nachbarnanteil,
                // für hohe Degree (r~1) dominiert der Degree selbst.
                buffer_score = (1 - r) * config.haa_theta * (static_cast<float>(cnt_adj_partitioned) / degree)
                              + r * (degree / static_cast<float>(config.d_max));
                break;
            }

            case BUFFER_SCORE_CMS: // Community - Majority Score
            {
                std::vector<int> hash_map(config.k, 0);
                int cnt_future_neighbors = 0;
                int most_common_partition_count = 0;
                for (LongNodeID adj_id : adjacents) {
                    PartitionID adj_part = (*config.stream_nodes_assign)[adj_id - 1];
                    if (adj_part != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                        hash_map[adj_part]++;
                        if (hash_map[adj_part] > most_common_partition_count) {
                            most_common_partition_count = hash_map[adj_part];
                            if (most_common_partition_count > degree / 2) {
                                break;
                            }
                        }
                    }
                }
                buffer_score =(float) most_common_partition_count /  degree; // +  (float) degree / config.d_max;
                break;
            }

            case BUFFER_SCORE_NSS: // Neighborhood Seen Score
            {
                unsigned cnt_adj_in_buffer = 0;
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    } else if (pq.contains(global_adj_id)) {
                    // } else if (node_id_to_buffer_item2.count(global_adj_id)) {
                        cnt_adj_in_buffer++;
                    }
                }
                // NOTE: Maybe account for degree somehow? higher importance to higher degree nodes exponentially or linearly? (e.g. degree / config.d_max)
                buffer_score = (float) (cnt_adj_partitioned + cnt_adj_in_buffer) /  degree;
                break;
            }
            case BUFFER_SCORE_GTS: // Graph Traversal Score
            {
                // NOTE: Maybe leave out this for performance reasons
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    }
                }
                buffer_score = (float) degree / config.d_max;
                break;
            }
            case BUFFER_SCORE_CBS3: // Cuttana buffer score 2
            {
                bool adj_is_partitioned;
                unsigned cnt_adj_in_buffer = 0;
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;

                    } else if (pq.contains(global_adj_id)) {
                    // } else if (node_id_to_buffer_item2.count(global_adj_id)) {
                        cnt_adj_in_buffer++;
                    }
                }
                // NOTE: on update: how to deal with cnt_adj_in_buffer?
                buffer_score = (float) degree / config.d_max + THETA * (float) (cnt_adj_partitioned + cnt_adj_in_buffer) /  degree; // Range: [0, 3] (first term: [0, 1], second term: [0, THETA])
                break;
            }
            case BUFFER_SCORE_CBS2: // Cuttana buffer score 2
            {
                bool adj_is_partitioned;
                unsigned cnt_adj_in_buffer = 0;
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    } else if (pq.contains(global_adj_id)) {

                    // } else if (node_id_to_buffer_item2.count(global_adj_id)) {
                        cnt_adj_in_buffer++;
                    }
                }
                // NOTE: on update: how to deal with cnt_adj_in_buffer?
                buffer_score = (float) degree / config.d_max + THETA * (float) (cnt_adj_partitioned + cnt_adj_in_buffer) /  degree; // Range: [0, 3] (first term: [0, 1], second term: [0, THETA])
                break;
            }
            case BUFFER_SCORE_CBS: // Cuttana buffer score
            default: // Default to classic buffer score
            {
                bool adj_is_partitioned;
                for (const LongNodeID& global_adj_id : adjacents) {
                    adj_is_partitioned = (*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION;
                    if (adj_is_partitioned) {
                        cnt_adj_partitioned++;
                    }
                }
                buffer_score = (float) degree /  config.d_max + THETA * (float) cnt_adj_partitioned / degree; // Range: [0, 3] (first term: [0, 1], second term: [0, THETA])
                break;
            }
        }


        return buffer_score;
    }

    float calc_updated_buffer_score(LongNodeID node_id, PQItem& buffer_item)  {
        const std::vector<LongNodeID>& adjacents = buffer_item.get_adjacents();
        float degree = static_cast<float>(adjacents.size());

        switch (config.buffer_score_type) {
            case BUFFER_SCORE_ANC:
                return buffer_item.buffer_score + 1.0 /  degree;
            case BUFFER_SCORE_ANC2:
                return buffer_item.buffer_score +  0.6 /  degree;
            case BUFFER_SCORE_HAA:
                {
                    float r = std::pow(degree / static_cast<float>(config.d_max), current_beta);
                    // float r = fast_pow(degree / config.d_max, current_beta);
                    // float r = fast_pow_specialized(degree / config.d_max, current_beta);
                    return std::min(buffer_item.buffer_score + (1 - r) * (1 / degree), 1.0f);
                }
            case BUFFER_SCORE_CMS:
                buffer_item.num_adj_partitioned = 0;
                // NOTE: only for evaluation, if it performs good, then some more sophisticated logic for updating should be implemented
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);
            case BUFFER_SCORE_NSS:
                buffer_item.num_adj_partitioned = 0;
                // NOTE: only for evaluation, if it performs good, then some more sophisticated logic for updating should be implemented
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);
            case BUFFER_SCORE_GTS:
                return MIN(buffer_item.buffer_score + config.gts_alpha, get_max_buffer_score(config));
            case BUFFER_SCORE_CBS2:
                buffer_item.num_adj_partitioned = 0;
                // NOTE: only for evaluation, if it performs good, then some more sophisticated logic for updating should be implemented
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);
            case BUFFER_SCORE_CBS:
            default:
                return buffer_item.buffer_score + (float) THETA / degree;
        }


    }



    float get_avg_degree() {
        return static_cast<float>(total_degree_sum) / node_counter;
    }

    void update_beta(double tau = 50.0) {
        float avg_degree = get_avg_degree();
        // current_beta = 1.5 * (1 + (avg_degree / config.d_max));
        switch (config.haa_hub_mode) {
            case HAA_INC_FOR_HUBS:
                current_beta = config.haa_beta_min + (config.haa_beta_max - config.haa_beta_min) * std::exp(-avg_degree / config.haa_tau);
                break;
            case HAA_DEC_FOR_HUBS:
                current_beta =  config.haa_beta_min - (config.haa_beta_max - config.haa_beta_min) * std::exp(-avg_degree / config.haa_tau);
                break;
        }
        // std::cout << "Updated beta: " << current_beta << std::endl;
    }

    // Adds a node to the buffer
    void addNode(LongNodeID global_node_id, std::vector<LongNodeID> &adjacents) {
        if (config.haa_hub_mode != HAA_NONADAPTIVE) {
            total_degree_sum += adjacents.size();
            node_counter++;
            if (node_counter == 1000) { // && node_counter % 100 == 0) {
                update_beta();
            }
        }

        // float buffer_score = calc_buffer_score(global_node_id, adjacents, node_id_to_buffer_item[global_node_id - 1].num_adj_partitioned);
        unsigned num_adj_partitioned = 0;
        float buffer_score = calc_buffer_score(global_node_id, adjacents, num_adj_partitioned);


        // node_id_to_buffer_item2[global_node_id] = PQItem(buffer_score, adjacents, num_adj_partitioned);
        // node_id_to_buffer_item2.emplace(
        //     global_node_id,
        //     PQItem(buffer_score, adjacents, num_adj_partitioned)
        // );
        // pq.insert(global_node_id,  discretize_score(buffer_score));

        // PQItem item(buffer_score, adjacents, num_adj_partitioned, 0, discretize_score(buffer_score));

        pq.insert(global_node_id, buffer_score, adjacents, num_adj_partitioned);
        // pq.insert(global_node_id, std::move(item));


        // node_id_to_buffer_item[global_node_id - 1].adjacents = new std::vector<LongNodeID>(adjacents);
        // node_id_to_buffer_item[global_node_id - 1].buffer_score = buffer_score;
    }

    // Removes and partitions the node with the highest priority
    void partitionTopNode() {
        LongNodeID node_id = pq.deleteMax();
        auto &buffer_item = pq.getBufferItem(node_id);

        assert((*config.stream_nodes_assign)[node_id - 1] == INVALID_PARTITION);

        // Partition the node
        // auto &adjacents = *node_id_to_buffer_item[node_id - 1].adjacents;
        // auto &adjacents = node_id_to_buffer_item2[node_id].get_adjacents();
        auto &adjacents = buffer_item.get_adjacents();
        // pq.getKey(node_id);
        partition_single_node(config, node_id, adjacents);

        // Update neighbors and clear buffer item
        update_neighbours_priority(adjacents);

        // node_id_to_buffer_item[node_id - 1].make_invalid();
        // node_id_to_buffer_item[node_id - 1].clean();
        // node_id_to_buffer_item2[node_id].make_invalid();
        completely_remove_node(node_id);
    }

    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesToBatch(std::vector<LongNodeID> *&input_idxs, LongNodeID batch_size) {
        // Initialize the partition configuration
        config.nmbNodes = MIN(batch_size, pq.size());
        input_idxs = new std::vector<LongNodeID>(config.nmbNodes);

        config.node_in_current_block_set = new google::dense_hash_set<LongNodeID>();
        config.node_in_current_block_set->set_empty_key(UNDEFINED_LONGNODE);
        config.node_in_current_block_set->resize(config.nmbNodes);

        // Extract the top batch_size number of nodes from the queue
        int node_counter = 0;
        while (node_counter < config.nmbNodes && !pq.empty()) {
            LongNodeID node_id = pq.deleteMax();
            assert((*config.stream_nodes_assign)[node_id - 1] == INVALID_PARTITION);

            auto &buffer_item = pq.getBufferItem(node_id);
            auto &adjacents = buffer_item.get_adjacents();
            // auto &buffer_item = node_id_to_buffer_item2[node_id];
            // pq.getKey(node_id);
            //  auto &buffer_item = node_id_to_buffer_item[node_id - 1];

            // Partition the node directly if it has a high degree
            if (adjacents.size() > config.d_direct) {
                partition_single_node(config, node_id, adjacents);

                // Update neighbors and make buffer item invalid
                update_neighbours_priority(adjacents, false);
                completely_remove_node(node_id);
                continue;
            }

            // Update neighbors and make buffer item invalid
            update_neighbours_priority(adjacents, false);

            // Add the node to the batch
            (*input_idxs)[node_counter] = node_id;
            config.node_in_current_block_set->insert(node_id);
            node_counter++;
        }

        if (node_counter < config.nmbNodes) {
            config.nmbNodes = node_counter;
            input_idxs->resize(node_counter);

        }

        if (config.haa_hub_mode != HAA_NONADAPTIVE) {
            update_beta();
        }
    }



    // Update the priority value of the neighbours of the node that was just partitioned in the priority queue
    void update_neighbours_priority(std::vector<LongNodeID> &adjacents, bool part_adj_directly = PARTITION_ADJ_DIRECTLY_ENABLED) {
        // Ensure that the adjacents contains at least one neighbor.
        if (adjacents.size() == 0) {
            return;
        }

        for (LongNodeID adj_id : adjacents) {
            bool is_partitioned = (*config.stream_nodes_assign)[adj_id - 1] != INVALID_PARTITION;

            // if (!is_partitioned && node_id_to_buffer_item[adj_id - 1].is_valid()) { // pq.contains(adj_id)) { // {
            if (!is_partitioned) {
                if (pq.contains_in_bucket(adj_id)) {
                // auto it = node_id_to_buffer_item2.find(adj_id);
                // if ( it != node_id_to_buffer_item2.end() && it->second.is_valid()) {
                    // auto &adj_buffer_item = it->second; // node_id_to_buffer_item2[adj_id];

                    PQItem& adj_buffer_item = pq.getBufferItem(adj_id);
                    // if (!is_partitioned && node_id_to_buffer_item2.count(adj_id) && node_id_to_buffer_item2[adj_id].is_valid()) {
                    // auto &adj_buffer_item = node_id_to_buffer_item[adj_id - 1];
                    // pq.getKey(adj_id);
                    auto &adj_adjacents = adj_buffer_item.get_adjacents();
                    int adj_degree = adj_adjacents.size();
                    if (adj_degree == 0) {
                        continue;
                    }

                    adj_buffer_item.num_adj_partitioned++;

                    // Check if all neighbours of the neighbour are partitioned, if so, partition the neighbour
                    if (part_adj_directly && adj_degree > 3 && adj_degree == adj_buffer_item.num_adj_partitioned) {
                        // assert_neighbors_partitioned(config, adj_buffer_item.adjacents, true);
                        // cnt_part_adj_directly++;
                        pq.deleteNode(adj_id);
                        partition_single_node(config, adj_id, adj_adjacents);

                        // Update neighbors and clear buffer item
                        update_neighbours_priority(adj_adjacents);

                        completely_remove_node(adj_id);
                    } else {
                        // Update buffer score of neighbours
                        // assert_neighbors_partitioned(config, adj_buffer_item.adjacents, false);
                        float updated_buffer_score = calc_updated_buffer_score(adj_id, adj_buffer_item);
                        // float score = discretize_score(updated_buffer_score);
                        pq.increaseKey(adj_id, updated_buffer_score);
                        // adj_buffer_item.buffer_score = updated_buffer_score;
                    }
                }
            }
        }
    }

    // Helper-Methoden
    bool isEmpty() const { return pq.empty(); }
    size_t size() const { return pq.size(); }

    std::vector<LongNodeID> &get_adjacents(LongNodeID node_id) {
        // return *node_id_to_buffer_item[node_id - 1].adjacents;
        // pq.getKey(node_id);

        return pq.getBufferItem(node_id).get_adjacents();

        // return node_id_to_buffer_item2[node_id].get_adjacents();

    }

    void completely_remove_node(LongNodeID node_id) {
        pq.completely_remove_node(node_id);
    }
};

#endif // BUFFER_NKJSAF9