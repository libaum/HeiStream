/******************************************************************************
 * init_fennel.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITFENNEL_7I4IR31Y
#define INITFENNEL_7I4IR31Y

#include "initial_partitioner.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "uncoarsening/refinement/refinement.h"
#include "definitions.h"
#include <algorithm>

struct pq_node {
    int id;
    float priority;

    pq_node(int id, float priority) : id(id), priority(priority) {}
};

struct ComparePriority {
    bool operator()(pq_node const &n1, pq_node const &n2) {
        // Higher priority first
        return n1.priority < n2.priority;
    }
};

class init_fennel : public initial_partitioner {
        public:
                init_fennel();
                virtual ~init_fennel();

                void initial_partition( PartitionConfig & config, const unsigned int seed,
                                graph_access & G, int* partition_map, int ismultisec=0);

                void initial_partition( PartitionConfig & config, const unsigned int seed,
                                graph_access & G,
                                int* xadj,
                                int* adjncy,
                                int* vwgt,
                                int* adjwgt,
                                int* partition_map,
                                int ismultisec=0);

        private:
		EdgeWeight fennel(PartitionConfig & partition_config, graph_access & G);

        float get_priority(graph_access &G, NodeID node, std::vector<bool> &node_is_partitioned, PartitionConfig &partition_config);
        void partition_node(PartitionConfig &partition_config, graph_access &G, NodeID node, std::vector<PartitionID> &hash_map, std::vector<NodeWeight> &cluster_sizes, std::vector<NodeWeight> &cluster_ghost_nodes, random_functions::fastRandBool<uint64_t> &random_obj, double fennel_weight, bool preliminary_sol, std::vector<bool> &node_is_partitioned);
        void update_neighbours_priority(graph_access &G, NodeID node_id, std::priority_queue<pq_node, std::vector<pq_node>, ComparePriority> &pq, std::vector<float> &node_id_to_priority, std::vector<bool> &node_is_partitioned, PartitionConfig &partition_config);
                                    //    random_functions::fastRandBool<uint64_t> &random_obj, double fennel_weight, bool preliminary_sol, std::vector<NodeWeight> &cluster_sizes, std::vector<NodeWeight> &cluster_ghost_nodes, std::vector<PartitionID> &hash_map);
};


#endif /* end of include guard: BIPARTITION_7I4IR31Y */
