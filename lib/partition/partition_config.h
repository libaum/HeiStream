/******************************************************************************
 * partition_config.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_CONFIG_DI1ES4T0
#define PARTITION_CONFIG_DI1ES4T0

#include "definitions.h"
#include "data_structure/buffered_map.h"
#include "data_structure/graph_access.h"
#include <map>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <unordered_set>

typedef struct {
        PartitionID block;
        int gain;
        int degree;
} DELTA;

class matrix;

// Configuration for the partitioning.
struct PartitionConfig
{
        PartitionConfig() {}

        LongNodeID max_block_weight;

        LongNodeID total_nodes_loaded;
        unsigned number_of_nodes;

        std::vector<NodeID> *local_to_global_map;
        std::vector<unsigned> *node_in_current_block;

        LongNodeID first_phase_buffer_len;
        LongNodeID second_phase_buffer_len;
        LongNodeID max_pq_size;
        int bq_disc_factor;

        //============================================================
        //=======================MATCHING=============================
        //============================================================
        bool edge_rating_tiebreaking;

        EdgeRating edge_rating;

        PermutationQuality permutation_quality;

        MatchingType matching_type;

        bool match_islands;

        bool first_level_random_matching;

        bool rate_first_level_inner_outer;

        NodeWeight max_vertex_weight;

        NodeWeight largest_graph_weight;

	NodeWeight work_load;

        unsigned aggressive_random_levels;

        bool disable_max_vertex_weight_constraint;

        //============================================================
        //===================INITIAL PARTITIONING=====================
        //============================================================
        unsigned int initial_partitioning_repetitions;

        unsigned int minipreps;

        bool refined_bubbling;

        InitialPartitioningType initial_partitioning_type;

        bool initial_partition_optimize;

        BipartitionAlgorithm bipartition_algorithm;

        bool initial_partitioning;

        int bipartition_tries;

        int bipartition_post_fm_limits;

        int bipartition_post_ml_limits;

        //============================================================
        //====================REFINEMENT PARAMETERS===================
        //============================================================
        bool corner_refinement_enabled;

        bool use_bucket_queues;

        RefinementType refinement_type;

        PermutationQuality permutation_during_refinement;

        ImbalanceType imbalance;

        unsigned bubbling_iterations;

        unsigned kway_rounds;

        bool quotient_graph_refinement_disabled;

        KWayStopRule kway_stop_rule;

        double kway_adaptive_limits_alpha;

        double kway_adaptive_limits_beta;

        unsigned max_flow_iterations;

        unsigned local_multitry_rounds;

        unsigned local_multitry_fm_alpha;

        bool graph_allready_partitioned;

        unsigned int fm_search_limit;

        unsigned int kway_fm_search_limit;

        NodeWeight upper_bound_partition;

        double bank_account_factor;

        RefinementSchedulingAlgorithm refinement_scheduling_algorithm;

        bool most_balanced_minimum_cuts;

        bool most_balanced_minimum_cuts_node_sep;

        unsigned toposort_iterations;

        bool softrebalance;

        bool rebalance;

        double flow_region_factor;

        bool gpa_grow_paths_between_blocks;

        //=======================================
        //==========GLOBAL SEARCH PARAMETERS=====
        //=======================================
        unsigned global_cycle_iterations;

        bool use_wcycles;

        bool use_fullmultigrid;

        unsigned level_split;

        bool no_new_initial_partitioning;

        bool omit_given_partitioning;

        StopRule stop_rule;

        int num_vert_stop_factor;

        bool no_change_convergence;

        //=======================================
        //===PERFECTLY BALANCED PARTITIONING ====
        //=======================================
	bool remove_negative_cycles;

        bool kaba_include_removal_of_paths;

        bool kaba_enable_zero_weight_cycles;

        double kabaE_internal_bal;

        CycleRefinementAlgorithm cycle_refinement_algorithm;

        int kaba_internal_no_aug_steps_aug;

        unsigned kaba_packing_iterations;

        bool kaba_flip_packings;

        MLSRule kaba_lsearch_p; // more localized search pseudo directed

        bool kaffpa_perfectly_balanced_refinement;

        unsigned kaba_unsucc_iterations;


        //=======================================
        //============PAR_PSEUDOMH / MH =========
        //=======================================
	double time_limit;

        double epsilon;

	unsigned no_unsuc_reps;

	unsigned local_partitioning_repetitions;

        bool mh_plain_repetitions;

        bool mh_easy_construction;

        bool mh_enable_gal_combine;

        bool mh_no_mh;

        bool mh_print_log;

        int  mh_flip_coin;

        int  mh_initial_population_fraction;

        bool mh_disable_cross_combine;

        bool mh_cross_combine_original_k;

        bool mh_disable_nc_combine;

        bool mh_disable_combine;

        bool mh_enable_quickstart;

        bool mh_disable_diversify_islands;

        bool mh_diversify;

        bool mh_diversify_best;

        bool mh_enable_tournament_selection;

        bool mh_optimize_communication_volume;

        unsigned mh_num_ncs_to_compute;

        unsigned mh_pool_size;

        bool combine; // in this case the second index is filled and edges between both partitions are not contracted

        unsigned initial_partition_optimize_fm_limits;

        unsigned initial_partition_optimize_multitry_fm_alpha;

        unsigned initial_partition_optimize_multitry_rounds;

        unsigned walshaw_mh_repetitions;

        unsigned scaleing_factor;

        bool scale_back;

	bool suppress_partitioner_output;

        unsigned maxT;

        unsigned maxIter;
        //=======================================
        //===============BUFFOON=================
        //=======================================
        bool disable_hard_rebalance;

        bool buffoon;

        bool kabapE;

        bool mh_penalty_for_unconnected;
        //=======================================
        //===============MISC====================
        //=======================================
        std::string input_partition;

        int seed;

        bool fast;

        bool eco;

        bool strong;

        bool kaffpaE;

	bool balance_edges;

        // number of blocks the graph should be partitioned in
        PartitionID k;

        bool compute_vertex_separator;

        bool only_first_level;

        bool use_balance_singletons;

        int amg_iterations;

        std::string graph_filename;

        std::string filename_output;

        bool kaffpa_perfectly_balance;

        bool mode_node_separators;

        //=======================================
        //===========SNW PARTITIONING============
        //=======================================
        NodeOrderingType node_ordering;

        int cluster_coarsening_factor;

        bool ensemble_clusterings;

        int label_iterations;

        int label_iterations_refinement;

        int number_of_clusterings;

        bool label_propagation_refinement;

        double balance_factor;

        bool cluster_coarsening_during_ip;

        bool set_upperbound;

        int repetitions;

        //=======================================
        //===========NODE SEPARATOR==============
        //=======================================
        int max_flow_improv_steps;

        int max_initial_ns_tries;

        double region_factor_node_separators;

	bool sep_flows_disabled;

	bool sep_fm_disabled;

	bool sep_loc_fm_disabled;

        int sep_loc_fm_no_snodes;

	bool sep_greedy_disabled;

	int sep_fm_unsucc_steps;

	int sep_loc_fm_unsucc_steps;

	int sep_num_fm_reps;

	int sep_num_loc_fm_reps;

        int sep_num_vert_stop;

        bool sep_full_boundary_ip;

        bool faster_ns;

        EdgeRating sep_edge_rating_during_ip;

        //=======================================
        //=========LABEL PROPAGATION=============
        //=======================================
        NodeWeight cluster_upperbound;

        //=======================================
        //=========INITIAL PARTITIONING==========
        //=======================================

        // variables controling the size of the blocks during
        // multilevel recursive bisection
        // (for the case where k is not a power of 2)
        std::vector<int> target_weights;

        bool initial_bipartitioning;

        int grow_target;

        //=======================================
        //===============QAP=====================
        //=======================================

        int communication_neighborhood_dist;

        LsNeighborhoodType ls_neighborhood;

        ConstructionAlgorithm construction_algorithm;

        DistanceConstructionAlgorithm distance_construction_algorithm;

        std::vector< int > group_sizes;

        std::vector< int > distances;

	int search_space_s;

        PreConfigMapping preconfiguration_mapping;

        int max_recursion_levels_construction;

        bool enable_mapping;


        //=======================================
        //===========integrated_mapping==========
        //=======================================

        bool integrated_mapping;
        bool multisection;
        bool qap_label_propagation_refinement;
        bool qap_blabel_propagation_refinement;
        bool qap_alabel_propagation_refinement;
        bool qap_multitry_kway_fm;
        bool qap_bmultitry_kway_fm;
        bool qap_kway_fm;
        bool qap_bkway_fm;
        bool qap_quotient_ref;
        bool qap_bquotient_ref;
        bool qap_0quotient_ref;
        bool bipartition_gp_local_search;
        bool skip_map_ls;
        bool suppress_output;
        bool no_change_convergence_map;
        bool full_matrix;
        matrix* D;
        std::vector< NodeID >* perm_rank;


        //=======================================
        //============= Delta gains =============
        //=======================================


        std::vector<std::pair<int,std::vector<DELTA*>>> *delta;
        std::vector<bool> *has_gains;
        bool use_delta_gains;
        bool quotient_more_mem;
        int *ref_layer;
        bool skip_delta_gains;



        //=======================================
        //======= Binary Online Distance ========
        //=======================================

        std::vector<std::vector<int>>  *bin_id;
        bool use_bin_id;
        std::vector<unsigned int>  *compact_bin_id;
        bool use_compact_bin_id;
        int bit_sec_len;
        int label_iterations_refinement_map;


        //=======================================
        //======== Adaptative Balancing =========
        //=======================================

        bool adapt_bal;
        double glob_block_upperbound;
        std::vector<int> interval_sizes;


        //=======================================
        //========== Stream Partition ===========
        //=======================================

        bool stream_input;
        LongNodeID stream_buffer_len;
        LongNodeID remaining_stream_nodes;
        LongEdgeID remaining_stream_edges;
        LongNodeID total_nodes;
        LongEdgeID total_edges;
        bool write_log;
        int remaining_stream_ew;
        LongNodeID total_stream_nodeweight;
        LongNodeID total_stream_nodecounter;
        LongNodeID stream_assigned_nodes;
        LongNodeID stream_n_nodes;
        std::ifstream* stream_in;
        LongNodeID lower_global_node;
        LongNodeID upper_global_node;
        std::vector<PartitionID>* stream_nodes_assign;
        std::vector<NodeWeight>* stream_blocks_weight;
        LongNodeID nmbNodes;
        std::vector<std::vector<EdgeWeight>> *degree_nodeBlock;
	/* std::vector<std::vector<NodeID>> *edge_block_nodes; */
	std::vector<std::vector<std::pair<NodeID,NodeWeight>>> *edge_block_nodes;
	int one_pass_algorithm;
        LongNodeID stream_total_upperbound;
        double fennel_gamma;
        double fennel_alpha;
        double fennel_alpha_gamma;
        bool use_fennel_objective; // maps global blocks to current stream blocks
        std::vector<NodeWeight>* add_blocks_weight;
	int fennel_dynamics;
        bool ram_stream;
        bool fennel_contraction;
	int fennel_batch_order;
	int quotient_nodes;
	int lhs_nodes;
        bool stream_initial_bisections;
	int n_batches;
	int curr_batch;
	double stream_global_epsilon;
        bool stream_output_progress;
	double batch_inbalance;
        bool skip_outer_ls;
	bool use_fennel_edgecut_objectives;

	// Initial partition via growing multiple BFS trees
	bool initial_part_multi_bfs;
	int multibfs_tries;

	// Initial partitioning via Fennel on the coarsest level
        int initial_part_fennel_tries;

	// Ghost neighbors
	buffered_map *ghostglobal_to_ghostkey;
	std::vector<NodeID>* ghostkey_to_node;
	std::vector<std::vector<std::pair<NodeID,ShortEdgeWeight>>>* ghostkey_to_edges;
	bool stream_whole_adjacencies;
	bool stream_allow_ghostnodes;
	LongNodeID ghost_nodes;
	int ghost_nodes_procedure;
	LongNodeID ghost_nodes_threshold;
	bool double_non_ghost_edges;

	// Restreaming and partial restreaming
	int num_streams_passes;
	int restream_number;
	bool restream_vcycle;

	int xxx;
    //=======================================
    //======= Stream Edge Partition ========
    //=======================================
    LongNodeID total_stream_edges;
    std::ofstream *stream_out;
    bool edge_partition;
    bool benchmark;
    bool dynamic_alpha;
    bool batch_alpha;
    bool minimal_mode;
    bool light_evaluator;
    bool convert_direct;
    bool use_queue;
    bool async_mode;
    bool evaluate_mode;
    bool include_weights;
    int parallel_nodes;
    int past_subset_size;
    int quotient_edges_count;
    double tau;
    NodeID num_split_edges;
    LongEdgeID remaining_stream_nodes_OG;
    LongEdgeID remaining_stream_graph_nodes;
    LongEdgeID fennel_edges;
    LongNodeID lower_global_node_conv;
    LongNodeID lower_global_store;
    unsigned long long start_pos;
    LongNodeID upper_global_node_conv;
    NodeID incremental_edge_ID;
    NodeID prev_batch_edge_ID;
    NodeID last_edge_count;
    NodeID back_node_count;
    NodeID forward_node_count;
    std::vector<std::vector<NodeID>> *nodes_on_edge_conv;
    std::vector<google::dense_hash_set<PartitionID>> *blocks_on_node;
    std::vector<NodeID> *blocks_on_node_minimal;
    double reps;
    double read_graph_time;
    double finding_past_assignments_time;
    double assign_edge_ID_time;
    double graph_model_time;
    double initial_partition_time;
    double stream_output_time;

        //=======================================
        // Conversion of graphs between formats =
        //=======================================

	bool relabel_nodes;
	int graph_translation_specs;
	bool input_header_absent;

        //=======================================
        //========== SpMxV partitioning =========
        //=======================================

	int matrix_m;
	int matrix_n;
	int matrix_nnz;


        //=======================================
        //===============Shared Mem OMP==========
        //=======================================
        bool enable_omp;

        void LogDump(FILE *out) const {
        }
};


#endif /* end of include guard: PARTITION_CONFIG_DI1ES4T0 */
