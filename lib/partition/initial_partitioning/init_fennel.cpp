/******************************************************************************
 * init_fennel.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "init_fennel.h"

init_fennel::init_fennel() {
}

init_fennel::~init_fennel() {
}

// This function performs the initial partitioning using the Fennel algorithm.
// It takes the following parameters:
// - config: the partitioning configuration
// - seed: the random seed for Fennel
// - G: the graph to be partitioned
// - partition_map: an array to store the partition assignments for each node
// - ismultisec: an optional parameter indicating whether to use multisector partitioning (default is 0)
void init_fennel::initial_partition(PartitionConfig &config,
                                    const unsigned int seed,
                                    graph_access &G,
                                    int *partition_map,
                                    int ismultisec /*=0*/) {

    timer t;
    t.restart();
    unsigned iterations = 1;
    EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();

    for (unsigned i = 0; i < iterations; i++) {
        fennel(config, G);
        G.set_partition_count(config.quotient_nodes);

        quality_metrics qm;
        EdgeWeight curcut = 0;

        // Calculate the cut value based on the Fennel objective or edge cut
        if (config.use_fennel_objective) {
            curcut = qm.fennel_objective(config, G, config.fennel_gamma, config.fennel_alpha);
        } else {
            curcut = qm.edge_cut(G);
        }

        // Update the best cut and partition map if a better cut is found
        if (curcut < best_cut) {
            best_cut = curcut;
            forall_nodes(G, n) {
                partition_map[n] = G.getPartitionIndex(n);
            }
            endfor
        }
    }
    PRINT(std::cout << "init_fennel took " << t.elapsed() << std::endl;)
}

void init_fennel::initial_partition(PartitionConfig &config,
                                    const unsigned int seed,
                                    graph_access &G,
                                    int *xadj,
                                    int *adjncy,
                                    int *vwgt,
                                    int *adjwgt,
                                    int *partition_map,
                                    int ismultisec /*=0*/) {

    std::cout << "not implemented yet" << std::endl;
}


float init_fennel::get_priority(graph_access &G, NodeID node, std::vector<bool> &node_is_partitioned, PartitionConfig &partition_config) {
    int partitioned_edge_weight = 0;
    int total_edge_weight = 0;
    int num_partioned_neighbours = 0;
    int num_total_neighbours = 0;
    int ghost_edge_weight = 0;
    int num_ghost_neighbours = 0;
    forall_out_edges(G, e, node) {
        num_total_neighbours++;
        EdgeWeight edge_weight = G.getEdgeWeight(e);
        total_edge_weight += edge_weight;

        NodeID target = G.getEdgeTarget(e);
        if (node_is_partitioned[target]) {
            partitioned_edge_weight += edge_weight;
            num_partioned_neighbours++;
        }
        if(target >= G.number_of_nodes() - partition_config.quotient_nodes) {
            ghost_edge_weight += edge_weight;
            num_ghost_neighbours++;
        }
    } endfor
    float theta = 2;
    float ghost_weight_factor = 1.0;

    // CBS
    return total_edge_weight / 300.0
        + theta * (partitioned_edge_weight / (float) total_edge_weight);

    // // CBS2 with total_edge_weight-quotient_edge_weight > 100
    // return total_edge_weight - ghost_edge_weight / 100.0
    //     + theta * ((partitioned_edge_weight-ghost_edge_weight) / (float) (total_edge_weight-ghost_edge_weight))
    //     + ghost_weight_factor * ghost_edge_weight / (float) partitioned_edge_weight;

    // CBS3
    // return total_edge_weight - ghost_edge_weight / 100.0
    //     + theta * ((partitioned_edge_weight-ghost_edge_weight) / (float) (total_edge_weight-ghost_edge_weight))
    //     + ghost_weight_factor * ghost_edge_weight / (float) partitioned_edge_weight
    //     + num_partioned_neighbours / (float) num_total_neighbours;

    // return total_edge_weight / 1000.0 + theta * (partitioned_edge_weight / (float) total_edge_weight + num_partioned_neighbours / (float) num_total_neighbours);
    // return total_edge_weight / 1000.0 + theta * partitioned_edge_weight / (float) total_edge_weight;
    // return (num_partioned_neighbours / (float) num_total_neighbours);// * (partitioned_edge_weight / (float) total_edge_weight);
}

// Update the priority value of the neighbours of the node that was just partitioned in the priority queue (insert with new value)
void init_fennel::update_neighbours_priority(graph_access &G, NodeID node_id, std::priority_queue<pq_node, std::vector<pq_node>, ComparePriority> &pq, std::vector<float> &node_id_to_priority, std::vector<bool> &node_is_partitioned, PartitionConfig &partition_config) {
    // random_functions::fastRandBool<uint64_t> &random_obj, double fennel_weight, bool preliminary_sol, std::vector<NodeWeight> &cluster_sizes, std::vector<NodeWeight> &cluster_ghost_nodes, std::vector<PartitionID> &hash_map) {
    forall_out_edges(G, e, node_id) {
        NodeID target = G.getEdgeTarget(e);
        if (!node_is_partitioned[target] && node_id_to_priority[target] != -1) {
            float new_priority = get_priority(G, target, node_is_partitioned, partition_config);
            // if ((int) new_priority == 1) {
            //     partition_node(partition_config, G, target, hash_map, cluster_sizes, cluster_ghost_nodes, random_obj, fennel_weight, preliminary_sol, node_is_partitioned);
            //     update_neighbours_priority(G, target, pq, node_id_to_priority, node_is_partitioned, partition_config,
            //         random_obj, fennel_weight, preliminary_sol, cluster_sizes, cluster_ghost_nodes, hash_map);
            //     continue;
            // }

            node_id_to_priority[target] = new_priority;
            pq.push(pq_node(target, new_priority));
        }
    } endfor
}


EdgeWeight init_fennel::fennel(PartitionConfig &partition_config, graph_access &G) {
    bool PARTITION_WITH_PQ = true;

    random_functions::fastRandBool<uint64_t> random_obj;
    // bool node_too_large = false;
    std::vector<PartitionID> hash_map(partition_config.k, 0);
    std::vector<NodeWeight> cluster_sizes(partition_config.k, 0);
    std::vector<NodeWeight> cluster_ghost_nodes(partition_config.k, 0);

    std::vector<bool> node_is_partitioned(G.number_of_nodes(), false);

    if (partition_config.restream_number) {
        for (NodeID node = 0; node < G.number_of_nodes() - partition_config.quotient_nodes; node++) {
            cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        }
    }
    for (NodeID node = G.number_of_nodes() - partition_config.quotient_nodes, end = G.number_of_nodes(); node < end; node++) {
        cluster_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);
        cluster_ghost_nodes[G.getPartitionIndex(node)] += G.getImplicitGhostNodes(node);
        node_is_partitioned[node] = true;
    }

    double fennel_weight = 2;
    double fennel_tmp = 0;
    switch (partition_config.fennel_dynamics) {
    case FENNELADP_ORIGINAL:
        fennel_weight = 1;
        break;
    case FENNELADP_DOUBLE:
        fennel_weight = 2;
        break;
    case FENNELADP_LINEAR:
        fennel_weight = 2 * partition_config.remaining_stream_nodes / (double)partition_config.stream_nodes_assign->size();
        break;
    case FENNELADP_MID_LINEAR:
        fennel_tmp = 2 * partition_config.remaining_stream_nodes / (double)partition_config.stream_nodes_assign->size();
        if (fennel_tmp <= 1) {
            fennel_weight = 2 * (fennel_tmp);
        }
        break;
    case FENNELADP_QUADRATIC:
        fennel_tmp = partition_config.remaining_stream_nodes / (double)partition_config.stream_nodes_assign->size();
        fennel_weight = 2 * fennel_tmp * fennel_tmp;
        break;
    case FENNELADP_MID_QUADRATIC:
        fennel_tmp = 2 * partition_config.remaining_stream_nodes / (double)partition_config.stream_nodes_assign->size();
        if (fennel_tmp <= 1) {
            fennel_weight = 2 * fennel_tmp * fennel_tmp;
        }
        break;
    case FENNELADP_MID_CONSTANT:
        fennel_tmp = partition_config.remaining_stream_nodes / (double)partition_config.stream_nodes_assign->size();
        if (fennel_tmp <= 1.5) {
            fennel_weight = 0.5;
        }
        break;
    case FENNELADP_EDGE_CUT:
        fennel_weight = 0;
        break;
    }

    //	if (partition_config.remaining_stream_nodes == 0 && partition_config.fennel_dynamics != FENNELADP_ORIGINAL) {
    //		fennel_weight = 0;
    //	}
    fennel_weight = 1;

    std::priority_queue<pq_node, std::vector<pq_node>, ComparePriority> pq_delayed_nodes;
    std::vector<float> node_id_to_priority(G.number_of_nodes(), -1); // if -1: not in priority queue, else: index in priority queue with corresponding priority

    for (int j = 0; j < partition_config.initial_part_fennel_tries; j++) {
        bool preliminary_sol = j || partition_config.restream_number;
		// Iterate over all nodes in the graph
		forall_nodes(G, node) {

            // Check if the node belongs to the quotient nodes (nodes that are not part of the remainder)
			if (node >= G.number_of_nodes() - partition_config.quotient_nodes) {
                continue;
			}

            if (PARTITION_WITH_PQ) {
                if (node_is_partitioned[node]) {
                    continue;
                }

                int num_of_neighbours = 0;
                int total_edge_weight = 0;
                int quotient_edge_weight = 0;
                forall_out_edges(G, e, node) {
                    num_of_neighbours++;
                    total_edge_weight += G.getEdgeWeight(e);
                    if (G.getEdgeTarget(e) >= G.number_of_nodes() - partition_config.quotient_nodes) {
                        quotient_edge_weight += G.getEdgeWeight(e);
                    }
                } endfor

                // First, assign nodes to the priority queue with priority based on the already partitioned neighbors or partition directly
                bool should_be_partitioned_directly = num_of_neighbours == 0 || total_edge_weight > 300; //total_edge_weight-quotient_edge_weight > 100;//  || total_edge_weight > 1000; //total_edge_weight-quotient_edge_weight > 70; //num_of_neighbours > 100 || (total_edge_weight > 200 && quotient_edge_weight/total_edge_weight > 0.7)
                if (should_be_partitioned_directly) {
                    // Partition node and update priority of neighbours in pq
                    partition_node(partition_config, G, node, hash_map, cluster_sizes, cluster_ghost_nodes, random_obj, fennel_weight, preliminary_sol, node_is_partitioned);
                    update_neighbours_priority(G, node, pq_delayed_nodes, node_id_to_priority, node_is_partitioned, partition_config);
                        // random_obj, fennel_weight, preliminary_sol, cluster_sizes, cluster_ghost_nodes, hash_map);
                } else if (node_id_to_priority[node] == -1) { // If node is not already in priority queue
                    // Add node to priority queue
                    float priority = get_priority(G, node, node_is_partitioned, partition_config);

                    // if (partitioned_edge_weight == total_edge_weight) {
                    //     partition_node(partition_config, G, node, hash_map, cluster_sizes, cluster_ghost_nodes, random_obj, fennel_weight, preliminary_sol, node_is_partitioned);
                    //     update_neighbours_priority(G, node, pq_delayed_nodes, node_id_to_priority, node_is_partitioned, partition_config,
                    //         random_obj, fennel_weight, preliminary_sol, cluster_sizes, cluster_ghost_nodes, hash_map);
                    //     continue;
                    // }

                    pq_delayed_nodes.push(pq_node(node, priority));
                    node_id_to_priority[node] = priority;
                }
            } else {
                partition_node(partition_config, G, node, hash_map, cluster_sizes, cluster_ghost_nodes, random_obj, fennel_weight, preliminary_sol, node_is_partitioned);
            }
        } endfor

        if (PARTITION_WITH_PQ) {
            while (!pq_delayed_nodes.empty()) {
                pq_node cur_node = pq_delayed_nodes.top();
                pq_delayed_nodes.pop();
                NodeID node = cur_node.id;
                if (!node_is_partitioned[node] && cur_node.priority == node_id_to_priority[node]) { // If priority is the current priority of the node
                    // Partition node and update priority of neighbours in pq

                    partition_node(partition_config, G, node, hash_map, cluster_sizes, cluster_ghost_nodes, random_obj, fennel_weight, preliminary_sol, node_is_partitioned);
                    update_neighbours_priority(G, node, pq_delayed_nodes, node_id_to_priority, node_is_partitioned, partition_config);
                        // random_obj, fennel_weight, preliminary_sol, cluster_sizes, cluster_ghost_nodes, hash_map);
                }
            }
        }

    }

    // for (bool node_id : node_is_partitioned) {
    //     if (!node_id) {
    //         std::cout << "Node not partitioned" << std::endl;
    //     }
    // }

    return 0;
}

void init_fennel::partition_node(PartitionConfig &partition_config, graph_access &G, NodeID node, std::vector<PartitionID> &hash_map, std::vector<NodeWeight> &cluster_sizes, std::vector<NodeWeight> &cluster_ghost_nodes, random_functions::fastRandBool<uint64_t> &random_obj, double fennel_weight, bool preliminary_sol, std::vector<bool> &node_is_partitioned) {
    bool node_too_large = true;
    double cur_value = 0;

    // Move the node to the cluster that is most common in the neighborhood
    forall_out_edges(G, e, node) {
        NodeID target = G.getEdgeTarget(e);
        // Exclude nodes that belong to the remainder or have already been processed in the preliminary solution
        // if (target >= G.number_of_nodes() - partition_config.quotient_nodes || target < node || preliminary_sol) {
        if (target >= G.number_of_nodes() - partition_config.quotient_nodes || node_is_partitioned[target] || preliminary_sol) {
            hash_map[G.getPartitionIndex(target)] += G.getEdgeWeight(e);
        }
    }
    endfor

    // Second sweep for finding the cluster with the maximum value and resetting the hash_map array
    PartitionID my_block = 0;
    if (preliminary_sol) {
        my_block = G.getPartitionIndex(node);
        cluster_sizes[my_block] -= G.getNodeWeight(node);
        cluster_ghost_nodes[my_block] -= G.getImplicitGhostNodes(node);
    }
    PartitionID max_block = my_block;
    double max_value = std::numeric_limits<double>::lowest();

    // Iterate over all clusters
    for (PartitionID cur_block = 0; cur_block < hash_map.size(); cur_block++) {
        cur_value = hash_map[cur_block];
        // Calculate the value for the current cluster
        cur_value -= fennel_weight * (G.getNodeWeight(node) * partition_config.fennel_alpha_gamma * random_functions::approx_sqrt(cluster_sizes[cur_block]));

        // Check if the current cluster has a higher value than the maximum value
        // If the values are equal, randomly choose one of them
        // Also check if the size of the cluster after adding the node is within the upperbound limit
        if ((cur_value > max_value || (cur_value == max_value && random_obj.nextBool())) && (cluster_sizes[cur_block] - cluster_ghost_nodes[cur_block] + G.getNodeWeight(node) - G.getImplicitGhostNodes(node) <= partition_config.stream_total_upperbound)) {
            node_too_large = false;
            max_value = cur_value;
            max_block = cur_block;
        }

        // Reset the value in the hash_map array
        hash_map[cur_block] = 0;
    }

    // If the node is too large for any cluster, assign it to a random cluster
    if (node_too_large) {
    max_block = fnv2a(partition_config.lower_global_node + node) % partition_config.k; // Random choice
    }

    // Update the cluster sizes and assign the node to the selected cluster
    cluster_sizes[max_block] += G.getNodeWeight(node);
    cluster_ghost_nodes[max_block] += G.getImplicitGhostNodes(node);
    G.setPartitionIndex(node, max_block);
    node_is_partitioned[node] = true;
}
