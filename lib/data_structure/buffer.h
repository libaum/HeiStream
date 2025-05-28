#ifndef BUFFER_NKJSAF9
#define BUFFER_NKJSAF9

#include <iostream>
#include <memory>
#include <optional>
#include <vector>
#include "timer.h"

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "macros_assertions.h"
#include "partition/partition_config.h"
#include "random_functions.h"

#define MIN(A, B) ((A) < (B)) ? (A) : (B)
#define MAX(A, B) ((A) > (B)) ? (A) : (B)

// CUTTANA HYPERPARAMETERS
const float THETA = 2;
// const int D_MAX = 1000;


inline void partition_single_node(PartitionConfig &partition_config, LongNodeID global_node_id, std::vector<LongNodeID> &adjacents) {

    partition_config.count_misc2++;
    std::vector<int> hash_map(partition_config.k, 0);
    int cnt_future_neighbors = 0;
    for (LongNodeID adj_id : adjacents) {
        PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id - 1];
        if (adj_part < partition_config.batch_manager->get_max_valid_partition_id()) {
            hash_map[adj_part]++;
        } else {
            cnt_future_neighbors++;
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

    // Assign the node to the partition with the best score
    (*partition_config.stream_nodes_assign)[global_node_id - 1] = best_partition;

    // Update partition load
    (*partition_config.stream_blocks_weight)[best_partition]++;


    if (partition_config.write_node_part_order) {
        if (!(*partition_config.node_part_order).empty() && (*partition_config.node_part_order).back().find(" -> ") == (*partition_config.node_part_order).back().size() - 4) {
            (*partition_config.node_part_order).back() += std::to_string(best_partition);
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
    double update_adj_time = 0.0;

public:
    Buffer(PartitionConfig &partition_config, LongNodeID max_pq_size)
        :   config(partition_config),
            pq(static_cast<unsigned>(std::floor(get_max_buffer_score(partition_config) * partition_config.bq_disc_factor)) + 1,
                partition_config.number_of_nodes, max_pq_size, partition_config.bq_disc_factor) {

        current_beta = config.haa_beta;
        total_degree_sum = 0;
        node_counter = 0;
        progress = 0.0;
    }


    static float get_max_buffer_score(const PartitionConfig& cfg) {
        switch (cfg.buffer_score_type) {
            case BUFFER_SCORE_ANR:
                return 1;
            case BUFFER_SCORE_HAA:
                return std::max(cfg.haa_theta, 1.0f);
            case BUFFER_SCORE_CBS:
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
            case BUFFER_SCORE_ANR: // Assigned Neighbor Count
            {
                for (const LongNodeID &global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
                        cnt_adj_partitioned++;
                    }
                }
                buffer_score = cnt_adj_partitioned /  degree;
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

            case BUFFER_SCORE_CBS: // Cuttana buffer score
            default: // Default to classic buffer score
            {
                for (const LongNodeID& global_adj_id : adjacents) {
                    if ((*config.stream_nodes_assign)[global_adj_id - 1] != INVALID_PARTITION) {
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
            case BUFFER_SCORE_ANR:
                return buffer_item.buffer_score + 1.0 /  degree;
            case BUFFER_SCORE_HAA:
                {
                    float r = std::pow(degree / static_cast<float>(config.d_max), current_beta);
                    // float r = fast_pow(degree / config.d_max, current_beta);
                    // float r = fast_pow_specialized(degree / config.d_max, current_beta);
                    return std::min(buffer_item.buffer_score + (1 - r) * (1 / degree), 1.0f);
                }
            case BUFFER_SCORE_CBS:
            default:
                return buffer_item.buffer_score + THETA / degree;
        }
    }

    float get_avg_degree() {
        return static_cast<float>(total_degree_sum) / node_counter;
    }

    void update_beta(double tau = 50.0) {
        // progress = static_cast<float>(node_counter) / config.number_of_nodes;
        // current_beta = 1.5 * (1 + (avg_degree / config.d_max));
        switch (config.haa_hub_mode) {
            case HAA_INC_FOR_HUBS:
                current_beta = config.haa_beta_min + (config.haa_beta_max - config.haa_beta_min) * std::exp(- get_avg_degree() / config.haa_tau);
                break;
            case HAA_DEC_FOR_HUBS:
                current_beta =  config.haa_beta_min - (config.haa_beta_max - config.haa_beta_min) * std::exp(- get_avg_degree() / config.haa_tau);
                break;
            case HAA_INC_FOR_PROG:
                current_beta = config.haa_beta_min +
                            progress * (config.haa_beta_max - config.haa_beta_min);
                break;

            case HAA_EXP_FOR_PROG:
                // Exponentiell ansteigend mit Fortschritt
                {
                    float exponent = config.haa_tau > 0.0 ? config.haa_tau : 2.0f;
                    current_beta = config.haa_beta_min +
                                    std::pow(progress, exponent) *
                                    (config.haa_beta_max - config.haa_beta_min);
                }
                break;

            case HAA_LOG_FOR_PROG:
                // Logarithmisch ansteigend (schnell am Anfang, langsam am Ende)
                {
                    // Logarithmische Kurve, die bei 0 beginnt und bei 1 endet
                    float log_prog = progress < 0.01f ? 0.0f :
                                    std::log(1.0f + 99.0f * progress) / std::log(100.0f);
                    current_beta = config.haa_beta_min +
                                  log_prog * (config.haa_beta_max - config.haa_beta_min);
                }
                break;

            case HAA_SIN_FOR_PROG:
                // Sinusförmig variierend (Welle)
                {
                    // Sinuswelle zwischen min und max, abhängig vom Fortschritt
                    // Komplett ein Zyklus über den gesamten Stream
                    float sin_val = (std::sin(progress * 2.0f * M_PI) + 1.0f) / 2.0f;
                    current_beta = config.haa_beta_min +
                                  sin_val * (config.haa_beta_max - config.haa_beta_min);
                }
                break;

            case HAA_SIGMOID_FOR_PROG:
                // S-förmig ansteigend (langsam am Anfang und Ende, schnell in der Mitte)
                {
                    // Sigmoid-Funktion, skaliert auf [0,1]
                    float x = (progress - 0.5f) * 10.0f; // Skalieren auf [-5, 5]
                    float sigmoid = 1.0f / (1.0f + std::exp(-x));
                    current_beta = config.haa_beta_min +
                                  sigmoid * (config.haa_beta_max - config.haa_beta_min);
                }
                break;

            case HAA_NONADAPTIVE:
                break;
        }
    }

    // Adds a node to the buffer or partition directly if buffer score is higher than max score
    bool addNode(LongNodeID global_node_id, std::vector<LongNodeID> &adjacents) {
        if (config.haa_hub_mode != HAA_NONADAPTIVE) {
            progress += 1.0 / config.number_of_nodes;
            node_counter++;
            total_degree_sum += adjacents.size();
            if (node_counter == 1000) { // && node_counter % 100 == 0) {
                update_beta();
            }
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

        if (config.write_node_part_order) {
            std::string entry = std::to_string(node_id) + " " + std::to_string(pq.getBufferItem(node_id).buffer_score) + " -> ";
            (*config.node_part_order).push_back(entry);
        }

        // Partition the node
        auto &adjacents = pq.getBufferItem(node_id).get_adjacents();
        partition_single_node(config, node_id, adjacents);

        // Update neighbors and clear buffer item
        update_neighbours_priority(adjacents);
        completely_remove_node(node_id);
    }

    bool max_value_below_cutoff() {
        // std::cout << "Max value: " << pq.maxValue() << " " << pq.maxValue() / 100.0f  << " " << config.bs_cutoff << std::endl;
        return pq.maxValue() / (float) config.bq_disc_factor < config.bs_cutoff;
    }

    bool max_value_above_1 (float val=0.0f) {
        return pq.maxValue() > val * config.bq_disc_factor;
    }

    // Loads the top nodes into a batch for MLP processing
    void loadTopNodesToBatch(std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>> *&batch_nodes, LongNodeID batch_size, size_t batch_id) {
        // Initialize the partition configuration

        batch_nodes = new std::vector<std::pair<LongNodeID, std::vector<LongNodeID>>>(MIN(batch_size, pq.size()));

        PartitionID batch_marker = config.batch_manager->get_batch_marker(batch_id);

        // LongNodeID min_batch_size = MIN(config.nmbNodes, 1000);
        // Extract the top batch_size number of nodes from the queue
        LongNodeID local_node_counter = 0;
        while (local_node_counter < config.nmbNodes && !pq.empty()) {
            // if (config.bs_cutoff != 0.0f) {
            //     if (local_node_counter > min_batch_size && max_value_below_cutoff()) {
            //         break;
            //     }
            // }

            LongNodeID node_id = pq.deleteMax();
            auto &adjacents = pq.getBufferItem(node_id).get_adjacents();

            if (config.write_node_part_order) {
                std::string entry = std::to_string(node_id) + " " + std::to_string(pq.getBufferItem(node_id).buffer_score) + " -> ";
                (*config.node_part_order).push_back(entry);
                // if (!(*partition_config.node_part_order).empty() && (*partition_config.node_part_order).back().find(" -> ") == (*partition_config.node_part_order).back().size() - 4) {
                //     (*partition_config.node_part_order).back() += std::to_string(best_partition);
                // }
            }

            // Partition the node directly if it has a high degree
            // if (adjacents.size() > config.d_direct) {
            //     partition_single_node(config, node_id, adjacents);

            //     // Update neighbors and make buffer item invalid
            //     update_neighbours_priority(adjacents, false);
            //     completely_remove_node(node_id);
            //     continue;
            // }

            (*config.stream_nodes_assign)[node_id - 1] = batch_marker;
            update_neighbours_priority(adjacents, false);



            (*batch_nodes)[local_node_counter] = std::make_pair(node_id, std::move(adjacents));
            completely_remove_node(node_id);

            local_node_counter++;



        }

        if (local_node_counter < batch_nodes->size()) {
            batch_nodes->resize(local_node_counter);
        }

        if (config.haa_hub_mode != HAA_NONADAPTIVE) {
            update_beta();
        }
    }



    // Update the priority value of the neighbours of the node that was just partitioned in the priority queue
    void update_neighbours_priority(std::vector<LongNodeID> &adjacents, bool part_adj_directly = false) {
        update_adj_t.restart();

        if (part_adj_directly == true) {
            part_adj_directly = config.part_adj_directly;
        }

        // Ensure that the adjacents contains at least one neighbor.
        if (adjacents.size() == 0) {
            return;
        }

        // Für Performance-Optimierung: direkter Zugriff
        auto& queue_map = pq.get_queue_index_map();
        std::unordered_map<LongNodeID, PQItem>::iterator it;

        for (LongNodeID adj_id : adjacents) {
            it = queue_map.find(adj_id);
            if (it != queue_map.end()) {
            // if (pq.contains_with_it(adj_id, it)) {
                PQItem& adj_buffer_item = it->second;
                auto &adj_adjacents = adj_buffer_item.get_adjacents();
                unsigned adj_degree = adj_adjacents.size();
                adj_buffer_item.num_adj_partitioned++;

                // Check if all neighbours of the neighbour are partitioned, if so, partition the neighbour
                if (part_adj_directly && adj_degree > 3 && adj_degree == adj_buffer_item.num_adj_partitioned ) { //&& config.buffer_score_type != BUFFER_SCORE_CBS2
                // if (part_adj_directly && adj_degree > config.param_int1 && adj_degree == adj_buffer_item.num_adj_partitioned && config.buffer_score_type != BUFFER_SCORE_CBS2) {
                    std::cout << "Partitioning directly: " << adj_id << " with degree: " << adj_degree << std::endl;
                    update_adj_time += update_adj_t.elapsed();
                    pq.deleteNode(adj_id);
                    partition_single_node(config, adj_id, adj_adjacents);

                    // Update neighbors and clear buffer item
                    update_neighbours_priority(adj_adjacents);
                    completely_remove_node(adj_id);
                    update_adj_t.restart();
                } else {
                    // Update buffer score of neighbours
                    float updated_buffer_score = calc_updated_buffer_score(adj_id, adj_buffer_item);
                    pq.increaseKey(adj_id, updated_buffer_score);
                }
            }
            // if ((*config.stream_nodes_assign)[adj_id - 1] == INVALID_PARTITION) {
            // }
        }
        update_adj_time += update_adj_t.elapsed();
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

    float get_max_value() {
        return pq.maxValue();
    }

    float deleteMax() {
        return pq.deleteMax();
    }

    void print_pq_statistics() const {
        pq.print_statistics();
    }
};

#endif // BUFFER_NKJSAF9