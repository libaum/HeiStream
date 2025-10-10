#ifndef BUFFER_NKJSAF9
#define BUFFER_NKJSAF9

#include <iostream>
#include <memory>
#include <optional>
#include <vector>
#include "timer.h"

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "macros_assertions.h"
#include "partition/partition_config.h"
#include "random_functions.h"
#include "batch_id_manager.h"

#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#define MAX(A, B) ((A) > (B)) ? (A) : (B)


inline void partition_single_node(PartitionConfig &partition_config, LongNodeID global_node_id, std::vector<LongNodeID> &adjacents) {

    partition_config.count_misc2++;
    std::vector<float> hash_map(partition_config.k, 0);
    for (LongNodeID adj_id : adjacents) {
        PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id - 1];
        if (adj_part < partition_config.k) {
        // if (adj_part < partition_config.batch_manager->get_max_valid_partition_id()) {
            hash_map[adj_part]++;
        } else if (adj_part < 2 * partition_config.k) {
        // if (adj_part < partition_config.batch_manager->get_max_valid_partition_id()) {
            hash_map[adj_part - partition_config.k]+= partition_config.ghost_importance; // Ghost nodes are only half as important
        }
    }

    // Iterate over partitions to compute FENNEL scores
    PartitionID best_partition = 0;
    float best_score = std::numeric_limits<float>::lowest();

    for (PartitionID cur_partition = 0; cur_partition < partition_config.k; ++cur_partition) {
        NodeWeight current_block_weight = (*partition_config.stream_blocks_weight)[cur_partition];
        // Skip or penalize partitions that are already "full"
        if (current_block_weight >= partition_config.max_block_weight) {
            // This partition is not feasible anymore
            continue;
        }

        float edge_gain = hash_map[cur_partition];
        NodeWeight partition_load = (*partition_config.stream_blocks_weight)[cur_partition];
        float load_penalty = partition_config.fennel_alpha_gamma * random_functions::approx_sqrt(partition_load);
        float score = edge_gain - load_penalty;

        if (score > best_score) {
            best_score = score;
            best_partition = cur_partition;
        }
    }

    // Assign the node to the partition with the best score
    (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;

    // Update partition load
    (*partition_config.stream_blocks_weight)[best_partition]++;

    if (partition_config.ghost_neighbors_enabled) {
        for (LongNodeID adj_id : adjacents) {
            PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id - 1];
            if ( adj_part == INVALID_PARTITION || partition_config.k <= adj_part && adj_part < 2 * partition_config.k ) {
                // If not yet assigned, assign temp partition ID of best_partition + k
                (*partition_config.stream_nodes_assign)[adj_id - 1] = best_partition + partition_config.k;
            }
        }
    }
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
    else if (exponent == 0.5f) {
        return random_functions::approx_sqrt(base);
    }
    else if (exponent == 0.0f) {
        return 1.0f;
    }
    // Allgemeiner Fall für andere Exponenten
    return std::exp(exponent * std::log(base));
}

class Buffer {
private:
    PartitionConfig &config;
    bucket_pq pq;

    LongNodeID total_degree_sum;
    LongNodeID node_counter;
    double progress;
    float current_beta;

    timer update_adj_t;
    TIMING_DECLARE(double update_adj_time) = 0.0;

    std::function<void(PartitionTask&&)> push_task_callback;

    float b_score_dmax; /// Maximum buffer score dmax, degree factor of all degrees higher than this is counted as 1.0 in the score calculation
    bool adaptive_theta_active;



public:
    Buffer(PartitionConfig &partition_config, LongNodeID max_pq_size)
        :   config(partition_config),
            pq(partition_config, static_cast<unsigned>(std::floor(get_max_buffer_score(partition_config) * partition_config.bq_disc_factor)) + 1,
                partition_config.number_of_nodes, max_pq_size, partition_config.bq_disc_factor) {

        b_score_dmax = std::min(static_cast<float>(config.d_max), 10000.0f); // Cap d_max to avoid too large values
        current_beta = config.haa_beta;

        adaptive_theta_active = config.haa_theta_min < config.haa_theta_max;
        total_degree_sum = 0;
        node_counter = 0;
        progress = 0.0;
    }

    void set_push_task_callback(std::function<void(PartitionTask&&)> callback) {
        push_task_callback = std::move(callback);
    }


    static float get_max_buffer_score(const PartitionConfig& cfg) {
        switch (cfg.buffer_score_type) {
            case BUFFER_SCORE_ANR:
            case BUFFER_SCORE_CMS:
            case BUFFER_SCORE_NSS:
            case BUFFER_SCORE_GTS:
                return 1;
            case BUFFER_SCORE_HAA:
                return std::max((std::max(1.0f, cfg.haa_theta)), cfg.haa_theta_max) + cfg.haa_theta0;
            case BUFFER_SCORE_HAA2:
            case BUFFER_SCORE_HAA3:
                return std::max(1.0f, cfg.haa_theta);
            case BUFFER_SCORE_CBS:
            case BUFFER_SCORE_CBSQ:
            default:
                return 1.0f + cfg.cbs_theta;
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


        /// Eventuell TODO: adjust all stream_nodes_assign to stream_nodes_batch_marker[global_adj_id - 1] != INVALID_PARTITION
        switch (config.buffer_score_type) {
            case BUFFER_SCORE_ANR: // Assigned Neighbor Count
            {
                for (const LongNodeID &global_adj_id : adjacents) {
                    PartitionID pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                    if (pid != INVALID_PARTITION) {
                        bool is_no_ghost = pid < config.k || pid > 2 * config.k;
                        if (is_no_ghost) {
                            cnt_adj_partitioned++;
                        }
                    }
                }
                buffer_score = cnt_adj_partitioned /  degree;
                break;
            }
            case BUFFER_SCORE_HAA:
            {
                PartitionID pid;
                if (!config.ghost_neighbors_enabled) {
                    for (const LongNodeID& global_adj_id : adjacents) {
                        pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                        if (pid != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        } else if (config.sep_batch_marker && (*config.stream_nodes_batch_marker)[global_adj_id - 1] != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        }
                    }
                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);

                    /// OPTIMIZATION for BETA = 2
                    // buffer_score = degree_factor * degree_factor
                    buffer_score = std::pow(degree_factor, current_beta)
                                    + (config.haa_theta0 + config.haa_theta * (1 - degree_factor)) * (static_cast<float>(cnt_adj_partitioned) / degree);

                    // buffer_score = (config.haa_theta0 + config.haa_theta * (1 - degree_factor)) * (static_cast<float>(cnt_adj_partitioned) / degree) + std::pow(degree_factor, current_beta); //std::pow(degree_factor, current_beta);
                } else {

                    unsigned cnt_adj_ghost = 0;
                    for (const LongNodeID& global_adj_id : adjacents) {
                        pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                        if (pid >= config.k && pid < 2 * config.k) {
                            cnt_adj_ghost++;
                        } else if (pid != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        } else if (config.sep_batch_marker && (*config.stream_nodes_batch_marker)[global_adj_id - 1] != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        }
                    }

                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);

                     /// OPTIMIZATION for BETA = 2
                    // buffer_score = degree_factor * degree_factor
                    buffer_score = std::pow(degree_factor, current_beta)
                                    + (config.haa_theta0 + config.haa_theta * (1 - degree_factor)) * (static_cast<float>(cnt_adj_partitioned) / degree
                                    + config.ghost_importance * (static_cast<float>(cnt_adj_ghost) / degree)); /// Add ghost neighbors to the score
                }

                break;
            }

            case BUFFER_SCORE_HAA2:
            {
                if (config.ghost_neighbors_enabled) {
                    unsigned cnt_adj_ghost = 0;
                    PartitionID pid;
                    for (const LongNodeID& global_adj_id : adjacents) {
                        pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                        if (pid >= config.k && pid < 2 * config.k) {
                            cnt_adj_ghost++;
                        } else if (pid != INVALID_PARTITION) {
                            // bool is_no_ghost = pid < config.k || pid > 2 * config.k;
                            // if (is_no_ghost) {
                                cnt_adj_partitioned++;
                            // } else {
                            //     cnt_adj_ghost++; // Count ghost neighbors
                            // }
                        } else if (config.sep_batch_marker && (*config.stream_nodes_batch_marker)[global_adj_id - 1] != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        }
                    }

                    /// Berechne den adaptiven Gewichtungsfaktor r(v)

                    // ALTERNATIVE IMPLEMENTATION
                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);
                    float r = std::pow(degree_factor, current_beta);

                    // OPTIMIZATION FOR BETA = 2
                    // float r = degree_factor * degree_factor; // Direct calculation avoiding intermediate variable


                    // float r = fast_pow(degree / b_score_dmax, current_beta);
                    // float r = fast_pow_specialized(degree / b_score_dmax, current_beta);

                    /// Berechne den Score: für niedrige Degree (r~0) dominiert der Nachbarnanteil, für hohe Degree (r~1) dominiert der Degree selbst.

                    // buffer_score = (1 - r) * ( config.haa_theta * (static_cast<float>(cnt_adj_partitioned) / degree)
                    //                         + config.ghost_importance * (static_cast<float>(cnt_adj_ghost) / degree)) /// Add ghost neighbors to the score
                    //                 + r * degree_factor;

                    buffer_score = degree_factor
                                    + config.haa_theta * (1 - degree_factor) * (static_cast<float>(cnt_adj_partitioned) / degree)
                                    + config.ghost_importance * (static_cast<float>(cnt_adj_ghost) / degree); /// Add ghost neighbors to the score

                } else {
                    PartitionID pid;
                    for (const LongNodeID& global_adj_id : adjacents) {
                        pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                        if (pid != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        } else if (config.sep_batch_marker && (*config.stream_nodes_batch_marker)[global_adj_id - 1] != INVALID_PARTITION) {
                            cnt_adj_partitioned++;
                        }
                    }
                    // Berechne den adaptiven Gewichtungsfaktor r(v)

                    /// Berechne den Score: für niedrige Degree (r~0) dominiert der Nachbarnanteil, für hohe Degree (r~1) dominiert der Degree selbst.
                    float degree_factor = std::min(degree / b_score_dmax, 1.0f);
                    float ANR = (static_cast<float>(cnt_adj_partitioned) / degree);

                    buffer_score = degree_factor + config.haa_theta * (1 - degree_factor) * ANR;
                }

                break;
            }


            case BUFFER_SCORE_HAA3:
            {
                PartitionID pid;
                for (const LongNodeID& global_adj_id : adjacents) {
                    pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                    if (pid != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    } else if (config.sep_batch_marker && (*config.stream_nodes_batch_marker)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    }
                }
                // x = D / float(D_MAX)
                // phi = 1 - ANR# **gamma
                // p = 2 - ANR
                // deg_term = (1.0 - phi)*x + phi*(x**p)

                // Score = deg_term + THETA1*(1.0 - deg_term)*ANR

                float deg_factor = std::min(degree / b_score_dmax, 1.0f);
                float ANR = (static_cast<float>(cnt_adj_partitioned) / degree);
                float phi = 1 - ANR;
                float p = 2 - ANR;
                float deg_term = (1.0f - phi) * deg_factor + phi * std::pow(deg_factor, p);
                buffer_score = deg_term + config.haa_theta * (1.0f - deg_term) * ANR;

                break;
            }



            case BUFFER_SCORE_CMS: // Community - Majority Score
            {
                std::vector<int> hash_map(config.k, 0);
                // int cnt_future_neighbors = 0;
                int most_common_partition_count = 0;
                for (LongNodeID adj_id : adjacents) {
                    PartitionID adj_part = (*config.stream_nodes_assign)[adj_id - 1];
                    if (adj_part < TO_BE_PARTITIONED - 10000) {
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
                buffer_score =(float) most_common_partition_count /  degree; // +  (float) degree / b_score_dmax;
                break;
            }

            case BUFFER_SCORE_NSS: // Neighborhood Seen Score
            {
                unsigned cnt_adj_in_buffer = 0;
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    } else if (pq.contains(global_adj_id)) {
                        cnt_adj_in_buffer++;
                    }
                }
                // NOTE: Maybe account for degree somehow? higher importance to higher degree nodes exponentially or linearly? (e.g. degree / b_score_dmax)
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
                buffer_score = std::min(degree / b_score_dmax, 1.0f);
                break;
            }

            case BUFFER_SCORE_CBS3: // Cuttana buffer score 3
            {
                unsigned cnt_adj_in_buffer = 0;
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    } else if (pq.contains(global_adj_id)) {

                        cnt_adj_in_buffer++;
                    }
                }
                buffer_score = std::min(degree / b_score_dmax, 1.0f) + config.cbs_theta * (float) cnt_adj_partitioned /  degree + (float) cnt_adj_in_buffer/ degree; // Range: [0, 3] (first term: [0, 1], second term: [0, config.cbs_theta])
                break;
            }

            case BUFFER_SCORE_CBSQ: // Cuttana buffer score 2
            {
                unsigned cnt_adj_ghost = 0;
                // Zähle, wie viele Nachbarn schon partitioniert sind.
                for (const LongNodeID& global_adj_id : adjacents) {
                    PartitionID pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                    if (pid != INVALID_PARTITION) {
                        bool is_no_ghost = pid < config.k || pid > 2 * config.k;
                        if (is_no_ghost) {
                            cnt_adj_partitioned++;
                        } else {
                            // Count ghost neighbors
                            cnt_adj_ghost++;
                        }
                    }
                }

                buffer_score = std::min(std::pow(degree / b_score_dmax, current_beta), 1.0f) + config.cbs_theta * (float) cnt_adj_partitioned / degree; // Range: [0, 3] (first term: [0, 1], second term: [0, config.cbs_theta])
                break;
            }

            case BUFFER_SCORE_CBS: // Cuttana buffer score
            {
                unsigned cnt_adj_ghost = 0;
                // Zähle, wie viele Nachbarn schon partitioniert sind.
                for (const LongNodeID& global_adj_id : adjacents) {
                    PartitionID pid = (*config.stream_nodes_assign)[global_adj_id - 1];
                    if (pid != INVALID_PARTITION) {
                        bool is_no_ghost = pid < config.k || pid > 2 * config.k;
                        if (is_no_ghost) {
                            cnt_adj_partitioned++;
                        } else {
                            // Count ghost neighbors
                            cnt_adj_ghost++;
                        }
                    }
                }

                buffer_score = std::min(degree / b_score_dmax, 1.0f) + config.cbs_theta * (float) cnt_adj_partitioned / degree; // Range: [0, 3] (first term: [0, 1], second term: [0, config.cbs_theta])
                // if (config.bscore_ghost) {
                //     // Add ghost neighbors to the score
                //     buffer_score += config.cbs_theta * config.ghost_importance * (float) cnt_adj_ghost / degree;
                // }
                break;
            }
            default: // Default to classic buffer score
                std::cerr << "Unknown buffer score type: " << config.buffer_score_type << std::endl;
                exit(1);

        }

        return buffer_score;
    }

    float calc_updated_buffer_score(LongNodeID node_id, PQItem& buffer_item, bool is_active_ghost_neighbor = false) {
        const std::vector<LongNodeID>& adjacents = buffer_item.get_adjacents();
        float degree = static_cast<float>(adjacents.size());

        switch (config.buffer_score_type) {
            case BUFFER_SCORE_ANR:
                return buffer_item.buffer_score + 1.0 /  degree;

            case BUFFER_SCORE_HAA:
                {
                    if (!config.ghost_neighbors_enabled) {
                        float degree_factor = std::min(degree / b_score_dmax, 1.0f);
                        return buffer_item.buffer_score + config.haa_theta * (1.0 - degree_factor) * (1.0 / degree);
                    } else {
                        // If neighboring node from which the update is triggered is an active ghost neighbor, adjust the score update accordingly
                        float ghost_neighbor_factor = is_active_ghost_neighbor ? (config.inv_ghost_importance) : 1.0f;
                        float degree_factor = std::min(degree / b_score_dmax, 1.0f);
                        return buffer_item.buffer_score + ghost_neighbor_factor * config.haa_theta * (1.0 - degree_factor) * (1.0 / degree);
                    }
                }

            case BUFFER_SCORE_HAA2:
                return buffer_item.buffer_score + config.haa_theta * (1.0 - std::min(degree / b_score_dmax, 1.0f)) * (1.0 / degree);

            case BUFFER_SCORE_HAA3:
                buffer_item.num_adj_partitioned = 0;
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);

            case BUFFER_SCORE_CMS:
                buffer_item.num_adj_partitioned = 0;
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);

            case BUFFER_SCORE_NSS:
                buffer_item.num_adj_partitioned = 0;
                return calc_buffer_score(node_id, adjacents, buffer_item.num_adj_partitioned);

            case BUFFER_SCORE_GTS:
                return MIN(buffer_item.buffer_score + config.gts_alpha, get_max_buffer_score(config));

            case BUFFER_SCORE_CBS:
            default:
                return buffer_item.buffer_score + config.cbs_theta / degree;
        }
    }

    float get_avg_degree() {
        return static_cast<float>(total_degree_sum) / node_counter;
    }

    // Adds a node to the buffer or partition directly if buffer score is higher than max score
    bool addNode(LongNodeID global_node_id, std::vector<LongNodeID> &adjacents) {

        if (adaptive_theta_active && node_counter % 1000 == 0) { // && node_counter % 100 == 0) {
            node_counter++;
            progress += 1000.0 / config.number_of_nodes;
            config.haa_theta = config.haa_theta_min + progress * (config.haa_theta_max - config.haa_theta_min);

        }


        unsigned num_adj_partitioned = 0;
        float buffer_score = calc_buffer_score(global_node_id, adjacents, num_adj_partitioned);
        unsigned bucket_idx = pq.discretize_score(buffer_score);

        // if (pq.size() >= config.max_pq_size && bucket_idx > pq.maxValue()) { // test1
        // if (pq.size() >= config.max_pq_size && bucket_idx >= pq.maxValue()) { // test2
        if (config.stream_buffer_len == 1 && pq.size() >= config.max_pq_size && bucket_idx >= pq.maxValue()) { // test3
            return false;
        } else {
            pq.insert(global_node_id, buffer_score, adjacents, num_adj_partitioned, bucket_idx);
            return true;
        }
    }

    // Removes and partitions the node with the highest priority
    void partitionTopNode() {
        LongNodeID node_id = pq.deleteMax();

        // Partition the node
        auto &adjacents = pq.getBufferItem(node_id).get_adjacents();
        partition_single_node(config, node_id, adjacents);

        // Update neighbors and clear buffer item
        update_neighbours_priority(adjacents);
        completely_remove_node(node_id);
    }

    bool max_value_above_1 (float val=0.0f) {
        return pq.maxValue() > val * config.bq_disc_factor;
    }

    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesToBatch(std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> *&batch_nodes, LongNodeID batch_size, size_t batch_id) {
        // Initialize the partition configuration

        config.nmbNodes = MIN(batch_size, pq.size());
        // batch_nodes = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>(config.nmbNodes);

        PartitionID batch_marker = config.batch_manager->get_batch_marker(batch_id);

        std::vector<PartitionID>& node_to_batch_marker = config.sep_batch_marker ? (*config.stream_nodes_batch_marker) : (*config.stream_nodes_assign);

        // Extract the top batch_size number of nodes from the queue
        LongNodeID local_node_counter = 0;
        while (local_node_counter < config.nmbNodes && !pq.empty()) {
            LongNodeID node_id = pq.deleteMax();
            auto &adjacents = pq.getBufferItem(node_id).get_adjacents();

            node_to_batch_marker[node_id - 1] = batch_marker;
            update_neighbours_priority(adjacents, false);

            batch_nodes->emplace_back(std::make_pair(node_id, std::move(adjacents)));
            completely_remove_node(node_id);

            // if (config.batch_extraction_strategy == BATCH_EXTRACTION_STRATEGY_COMPLETE_BATCH_WITH_ADJ) {
            //     for (LongNodeID &adj_id : adjacents) {
            //         if (pq.contains(adj_id) && local_node_counter < config.nmbNodes - 2) {
            //             pq.deleteNode(adj_id);
            //             std::vector<LongNodeID> adj_adjacents = std::move(get_adjacents(adj_id));
            //             (*config.stream_nodes_assign)[adj_id - 1] = batch_marker;

            //             update_neighbours_priority(adj_adjacents, false);
            //             completely_remove_node(adj_id);
            //             // batch_nodes.emplace_back(adj_id, std::move(adj_adjacents));
            //             (*batch_nodes)[local_node_counter] = std::make_pair(adj_id, std::move(adj_adjacents));

            //             local_node_counter++;
            //         }
            //     }
            // }

            local_node_counter++;
        }

    }


    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesAndNeighborsToBatch(std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> &batch_nodes, LongNodeID batch_size, size_t batch_id) {
        // Initialize the partition configuration

        // batch_nodes = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>(MIN(batch_size, pq.size()));
        // batch_nodes->clear();

        PartitionID batch_marker = config.batch_manager->get_batch_marker(batch_id);
        std::vector<PartitionID>& node_to_batch_marker = config.sep_batch_marker ? (*config.stream_nodes_batch_marker) : (*config.stream_nodes_assign);
        // LongNodeID min_batch_size = MIN(config.nmbNodes, 1000);
        // Extract the top batch_size number of nodes from the queue
        LongNodeID local_node_counter = 0;
        while (local_node_counter < config.nmbNodes && !pq.empty()) {

            LongNodeID node_id = pq.deleteMax();
            std::vector<LongNodeID> adjacents = std::move(get_adjacents(node_id));

            node_to_batch_marker[node_id - 1] = batch_marker;
            // (*config.stream_nodes_assign)[node_id - 1] = batch_marker;
            // (*config.stream_nodes_batch_marker)[node_id - 1] = batch_marker;

            update_neighbours_priority(adjacents, false);
            completely_remove_node(node_id);

            if (config.batch_extraction_strategy == BATCH_EXTRACTION_STRATEGY_COMPLETE_BATCH_WITH_ADJ) {
                for (LongNodeID &adj_id : adjacents) {
                    if (pq.contains(adj_id) && local_node_counter < config.nmbNodes - 2) {
                        pq.deleteNode(adj_id);
                        std::vector<LongNodeID> adj_adjacents = std::move(get_adjacents(adj_id));
                        node_to_batch_marker[adj_id - 1] = batch_marker;
                        // (*config.stream_nodes_assign)[adj_id - 1] = batch_marker;
                        // (*config.stream_nodes_batch_marker)[adj_id - 1] = batch_marker;

                        update_neighbours_priority(adj_adjacents, false);
                        completely_remove_node(adj_id);
                        batch_nodes.emplace_back(adj_id, std::move(adj_adjacents));

                        local_node_counter++;
                    }
                }
            }

            batch_nodes.emplace_back(node_id, std::move(adjacents));

            local_node_counter++;
        }

        if (local_node_counter < batch_nodes.size()) {
            batch_nodes.resize(local_node_counter);
        }

    }


    // Update the priority value of the neighbours of the node that was just partitioned in the priority queue
    void update_neighbours_priority(std::vector<LongNodeID> &adjacents, bool part_adj_directly = false) {
        TIMING_START(update_adj_t);

        if (part_adj_directly == true) {
            part_adj_directly = config.part_adj_directly;
        }

        // Ensure that the adjacents contains at least one neighbor.
        if (adjacents.size() == 0) {
            return;
        }

        // Für Performance-Optimierung: direkter Zugriff
        // auto& queue_map = pq.get_queue_index_map();
        // std::unordered_map<LongNodeID, PQItem>::iterator it;

        for (LongNodeID adj_id : adjacents) {
            // if ((*config.stream_nodes_assign)[adj_id - 1] == INVALID_PARTITION) { /// TODO: THIS IS A BUG I THINK, because a node can be assigned to a ghost partition
                if (pq.contains(adj_id)) {
                    PQItem& adj_buffer_item = pq.getBufferItem(adj_id);
                    auto &adj_adjacents = adj_buffer_item.get_adjacents();
                    unsigned adj_degree = adj_adjacents.size();
                    adj_buffer_item.num_adj_partitioned++;

                    // Check if all neighbours of the neighbour are partitioned, if so, partition the neighbour
                    if (part_adj_directly && adj_degree > 3 && adj_degree == adj_buffer_item.num_adj_partitioned ) {
			TIMING_ACCUMULATE(update_adj_time, update_adj_t);
                        pq.deleteNode(adj_id);
                        partition_single_node(config, adj_id, adj_adjacents);

                        // Update neighbors and clear buffer item
                        update_neighbours_priority(adj_adjacents);
                        completely_remove_node(adj_id);
		        TIMING_START(update_adj_t);
                    } else {
                        // Update buffer score of neighbours
                        float updated_buffer_score = calc_updated_buffer_score(adj_id, adj_buffer_item);
                        pq.increaseKey(adj_id, updated_buffer_score);
                    }
                }
            // }
        }
        TIMING_ACCUMULATE(update_adj_time, update_adj_t);
    }

    // Update the priority value of the neighbours of the node that was just partitioned in the priority queue
    void update_neighbours_priority_parallel(std::vector<LongNodeID> &adjacents,
                                             bool part_adj_directly = false,
                                             bool is_active_ghost_neighbor = false) {

        TIMING_START(update_adj_t);


        if (part_adj_directly == true) {
            part_adj_directly = config.part_adj_directly;
        }

        // Ensure that the adjacents contains at least one neighbor.
        if (adjacents.size() == 0) {
            return;
        }

        // Für Performance-Optimierung: direkter Zugriff
        // auto& queue_map = pq.get_queue_index_map();
        // std::unordered_map<LongNodeID, PQItem>::iterator it;

        for (LongNodeID adj_id : adjacents) {
            // it = queue_map.find(adj_id);
            // if (it != queue_map.end()) {
            // if (pq.contains_with_it(adj_id, it)) {
            if (pq.contains(adj_id)) {
                // PQItem& adj_buffer_item = it->second;
                PQItem& adj_buffer_item = pq.getBufferItem(adj_id); // Ensure the item is loaded into the buffer
                auto &adj_adjacents = adj_buffer_item.get_adjacents();
                unsigned adj_degree = adj_adjacents.size();
                adj_buffer_item.num_adj_partitioned++;

                // Check if all neighbours of the neighbour are partitioned, if so, partition the neighbour
                if (part_adj_directly && adj_degree > 3 && adj_degree == adj_buffer_item.num_adj_partitioned ) {

                    // move adj_adjacents to a new vector to avoid dangling references
                    std::vector<LongNodeID> adj_adjacents_copy = std::move(adj_adjacents);

                    pq.deleteNode(adj_id);
                    completely_remove_node(adj_id);

                    bool adj_is_active_ghost_neighbor = false;
                    if (config.ghost_neighbors_enabled) {
                        // Check if adj is active ghost neighbor
                        adj_is_active_ghost_neighbor = config.k <= (*config.stream_nodes_assign)[adj_id - 1]
                                                        && (*config.stream_nodes_assign)[adj_id - 1] < 2 * config.k;
                    }

                    update_neighbours_priority_parallel(
                        adj_adjacents_copy,
                        true,
                        adj_is_active_ghost_neighbor
                    );
                    (*config.stream_nodes_assign)[adj_id - 1] = TO_BE_PARTITIONED;

                    PartitionTask task(
                        -1,
                        std::vector<BatchNode>{{adj_id, std::move(adj_adjacents_copy)}}
                    );

                    if (push_task_callback) {
                        push_task_callback(std::move(task));
                    }

                } else {

                    // Update buffer score of neighbours
                    float updated_buffer_score = calc_updated_buffer_score(adj_id, adj_buffer_item, is_active_ghost_neighbor);
                    pq.increaseKey(adj_id, updated_buffer_score);
                }
            }
        }
        TIMING_ACCUMULATE(update_adj_time, update_adj_t);
    }


    double get_update_adj_time() {
        return update_adj_time;
    }

    // Helper-Methoden
    bool isEmpty() const { return pq.empty(); }
    size_t size() const { return pq.size(); }

    std::vector<LongNodeID> &get_adjacents(LongNodeID node_id) {
        return pq.getBufferItem(node_id).get_adjacents();
    }

    void completely_remove_node(LongNodeID node_id) {
        pq.completely_remove_node(node_id);
    }

    unsigned get_max_value() {
        return pq.maxValue();
    }

    LongNodeID deleteMax() {
        return pq.deleteMax();
    }

    void print_pq_statistics() const {
        pq.print_statistics();
    }
};

#endif // BUFFER_NKJSAF9
