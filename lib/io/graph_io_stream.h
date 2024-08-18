/******************************************************************************
 * graph_io_stream.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPHIOSTREAM_H_
#define GRAPHIOSTREAM_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <list>
#include <algorithm>

#include "definitions.h"
#include "data_structure/graph_access.h"
#include "partition/partition_config.h"
#include "timer.h"
#include "random_functions.h"
#include "data_structure/buffered_map.h"

typedef std::vector<std::string>* LINE_BUFFER;

class graph_io_stream {
        public:
                graph_io_stream();
                virtual ~graph_io_stream () ;

                static
		NodeID createModel (PartitionConfig & config, graph_access & G, std::vector<std::vector<LongNodeID>>* &input);

                static
		void processNodeWeight(PartitionConfig & config, std::vector<NodeWeight>& all_nodes, NodeID node, NodeWeight weight, LongNodeID global_node);

                static
                void generalizeStreamPartition(PartitionConfig & config, graph_access & G_local);

                static
		void countAssignedNodes(PartitionConfig & config);

                static
                void onePassPartition(PartitionConfig & config, std::vector<std::vector<EdgeWeight>> & edges_virtualReal,
					std::vector<PartitionID> & blockVirtualToReal, std::vector<NodeWeight> & weight_VirtualBlocks);

		static
		int onePassDecide(PartitionConfig & config, NodeID node, std::vector<EdgeWeight> & edges_i_real);

		static
		double getFennelWeight(PartitionConfig & partition_config);

                static
		void writePartitionStream(PartitionConfig & config, const std::string & filename);

				static
 		float graph_io_stream::get_ratio_of_partitioned_neighbours(PartitionConfig & partition_config, std::vector<LongNodeID>* line, int lower_global_node, int max_global_node_in_batch, bool read_nw, bool read_ew);

                static
		void readFirstLineStream(PartitionConfig & partition_config, std::string graph_filename, EdgeWeight& total_edge_cut);

                static
		void loadRemainingLines(PartitionConfig & partition_config, LINE_BUFFER &lines);

                static
		void loadBufferLines(PartitionConfig & partition_config, LINE_BUFFER &lines, LongNodeID num_lines);

                static
		std::vector<std::string>* loadLinesFromStream(PartitionConfig & partition_config, LongNodeID num_lines);

                static
		void sortBatchByDegree(PartitionConfig & config);

                static
		void createGraphForBatch(PartitionConfig & config, graph_access & G, NodeID node_counter, EdgeID edge_counter,
			std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes);

                static
		void recoverBlockAssignedToNode(PartitionConfig & config, graph_access & G, NodeID node, NodeID node_counter);

                static
		void setupForGhostNeighbors(PartitionConfig & config);

                static
		void processGhostNeighborInBatch(PartitionConfig & config, NodeID node, LongNodeID ghost_target, EdgeWeight edge_weight);

                static
		void processQuotientEdgeInBatch(PartitionConfig & config, NodeID node, LongNodeID global_target, EdgeWeight edge_weight);

                static
		EdgeID insertRegularEdgeInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges,
										NodeID node, NodeID target, EdgeWeight edge_weight);

                static
		NodeID mapGhostKeysToNodesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges,
							std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes, NodeID& node_counter);

                static
		NodeID restreamMapGhostKeysToNodes(PartitionConfig & config);

                static
		NodeID greedyMapGhostKeysToNodes(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges,
						std::vector<NodeWeight>& all_nodes, std::vector<NodeWeight>& all_assigned_ghost_nodes, NodeID& node_counter);

                static
		EdgeID insertGhostEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edgesInBatch);

                static
		void insertQuotientNodesInBatch(PartitionConfig & config, std::vector<NodeWeight>& all_nodes, NodeID uncontracted_ghost_nodes, NodeID& node_counter);

                static
		EdgeID insertQuotientEdgesInBatch(PartitionConfig & config, std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID uncontracted_ghost_nodes);


                static
		EdgeID includeEdgeInBatch(std::vector<std::vector<std::pair<NodeID,EdgeWeight>>>& all_edges, NodeID node, NodeID target, EdgeWeight edge_weight);

                static
		void prescribeBufferInbalance(PartitionConfig & partition_config);

                static
		void streamEvaluatePartition(PartitionConfig & config, const std::string & filename, EdgeWeight& edgeCut);

                static
		void loadRemainingLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input, std::deque<std::vector<LongNodeID>> &delayed_lines_queue);

                static
		void loadBufferLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input, LongNodeID num_lines, std::deque<std::vector<LongNodeID>> &delayed_lines_queue);

                static
		std::vector<std::vector<LongNodeID>>* loadLinesFromStreamToBinary(PartitionConfig & partition_config, LongNodeID num_lines, std::deque<std::vector<LongNodeID>> &delayed_lines_queue);

		template< typename T>
                static
		T return_and_delete_element(std::vector<T> & vec, LongNodeID pos);

//                static
//		int readEdgeStream_writeMetisBuffered(const std::string & graph_filename, std::string filename_output, bool relabel_nodes);

};

inline void graph_io_stream::loadRemainingLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input, std::deque<std::vector<LongNodeID>> &delayed_lines_queue) {
	if (partition_config.ram_stream) {
		partition_config.max_delayed_nodes = 0; // no delayed nodes in RAM mode
		input = graph_io_stream::loadLinesFromStreamToBinary(partition_config, partition_config.remaining_stream_nodes, delayed_lines_queue);
	}
}

inline void graph_io_stream::loadBufferLinesToBinary(PartitionConfig & partition_config, std::vector<std::vector<LongNodeID>>* &input, LongNodeID num_lines, std::deque<std::vector<LongNodeID>> &delayed_lines_queue) {
	if (!partition_config.ram_stream) {
		input = graph_io_stream::loadLinesFromStreamToBinary(partition_config, num_lines, delayed_lines_queue);
	} else {
		std::fill(partition_config.node_in_current_block->begin(), partition_config.node_in_current_block->end(), 0);
		int lower_global_node = partition_config.number_of_nodes - partition_config.remaining_stream_nodes;
		int upper_global_node = lower_global_node + num_lines > partition_config.number_of_nodes ? partition_config.number_of_nodes : lower_global_node + num_lines;
		for (int global_node_id = lower_global_node+1; global_node_id <= upper_global_node; global_node_id++) {
			(*partition_config.node_in_current_block)[global_node_id-1] = 1;
		}
	}
}

inline float graph_io_stream::get_ratio_of_partitioned_neighbours(PartitionConfig & partition_config, std::vector<LongNodeID>* line, int lower_global_node, int max_global_node_in_batch, bool read_nw, bool read_ew) {
	// check if node should be delayed
	unsigned num_neighbours = 0;
	unsigned num_of_neighours_partitioned = 0;
	unsigned num_neighbours_in_batch = 0;
	unsigned col_counter = read_nw ? 2 : 1;
	while (col_counter < line->size()) {
		LongNodeID target = (*line)[col_counter++];
		if(read_ew) {
			col_counter++;
		}
		if (!partition_config.curr_batch == 0 && (*partition_config.stream_nodes_assign)[target-1] != INVALID_PARTITION) {
			num_of_neighours_partitioned++;
		}
		if ((*partition_config.node_in_current_block)[target-1] == 1 || (lower_global_node < target && target < max_global_node_in_batch)) {
			num_neighbours_in_batch++;
		}
		num_neighbours++;
	}
	if (num_neighbours == 0) {
		return 1;
	} else {
		return (num_of_neighours_partitioned + num_neighbours_in_batch) / (float) num_neighbours;
	}
}

inline std::vector<std::vector<LongNodeID>>* graph_io_stream::loadLinesFromStreamToBinary(PartitionConfig & partition_config, LongNodeID num_lines, std::deque<std::vector<LongNodeID>> &delayed_lines_queue) {
	bool read_ew = false;
	bool read_nw = false;
	std::vector<std::vector<LongNodeID>>* input;
	input = new std::vector<std::vector<LongNodeID>>(num_lines);
	std::vector<std::string>* lines;
	lines = new std::vector<std::string>(1);
	LongNodeID node_counter = 0;
	buffered_input *ss2 = NULL;
	std::fill(partition_config.node_in_current_block->begin(), partition_config.node_in_current_block->end(), 0);

	bool is_last_batch = partition_config.remaining_stream_nodes == partition_config.nmbNodes;
	bool is_second_last_batch = partition_config.remaining_stream_nodes <= 2*partition_config.nmbNodes;
	int max_capacity_delayed_nodes = is_second_last_batch ? (partition_config.remaining_stream_nodes - partition_config.nmbNodes < partition_config.max_delayed_nodes ? partition_config.remaining_stream_nodes - partition_config.nmbNodes : partition_config.max_delayed_nodes) : partition_config.max_delayed_nodes;

	int num_nodes_delayed = delayed_lines_queue.size();
	int lower_global_node = partition_config.total_nodes_loaded;
	int max_global_node_in_batch = partition_config.total_nodes_loaded + num_lines - node_counter;


	switch(partition_config.remaining_stream_ew) {
		case 1:
			read_ew = true;
			break;
		case 10:
			read_nw = true;
			break;
		case 11:
			read_ew = true;
			read_nw = true;
			break;
	}

	// Load new nodes until batch is full, delay nodes if criteria is met
	while( node_counter < num_lines) {
		if (num_nodes_delayed > 0) { // Before loading new nodes, check if delayed nodes can be processed
			std::vector<LongNodeID> *line = &delayed_lines_queue.front();
			// Check if line should stay delayed -> append to delayed nodes else append to input
			float ratio = get_ratio_of_partitioned_neighbours(partition_config, line, INFINITY, -INFINITY, read_nw, read_ew);
			bool should_stay_delayed = ratio <= partition_config.threshold_delay;
			if (should_stay_delayed) {
				delayed_lines_queue.push_back(*line);
				delayed_lines_queue.pop_front();
			} else {
				(*input)[node_counter++] = (*line);
				(*partition_config.node_in_current_block)[(*line)[0]-1] = 1; // line[0]: global_node_id
				delayed_lines_queue.pop_front();
			}
			num_nodes_delayed--;
			if (num_nodes_delayed == 0) {
				max_global_node_in_batch = partition_config.total_nodes_loaded + num_lines - node_counter;
			}
		} else if (!partition_config.stream_in->eof()) { // Load new nodes
			std::getline(*(partition_config.stream_in),(*lines)[0]);
			if ((*lines)[0][0] == '%') { // a comment in the file
				continue;
			}
			ss2 = new buffered_input(lines);
			LongNodeID global_node_id = partition_config.total_nodes_loaded+1;
			std::vector<LongNodeID> *new_line = &(*input)[node_counter];
			new_line->clear();
			new_line->push_back(global_node_id);
			ss2->simple_scan_line((*new_line), false);

			// Check if line should be delayed -> append to delayed nodes else append to input
			bool should_be_delayed = false;
			if (!is_last_batch) {
				bool delayed_nodes_has_capacity = delayed_lines_queue.size() < max_capacity_delayed_nodes;
				if (delayed_nodes_has_capacity) {
					float ratio = get_ratio_of_partitioned_neighbours(partition_config, new_line, lower_global_node, max_global_node_in_batch, read_nw, read_ew);
					should_be_delayed = ratio <= partition_config.threshold_delay;
				}
			}

			if (should_be_delayed) {
				delayed_lines_queue.push_back(*new_line);
			} else {
				if (!partition_config.ram_stream) {
					(*partition_config.node_in_current_block)[global_node_id-1] = 1;
				}
				node_counter++;
			}

			(*lines)[0].clear(); delete ss2;
			partition_config.total_nodes_loaded++;
		} else { // If the file has ended, we need to check if there are still delayed nodes to be processed
			if (delayed_lines_queue.size() == 0) {
				break;
			}
			// No check if line should be delayed -> just append to input
			std::vector<LongNodeID> *line = &delayed_lines_queue.front();
			(*input)[node_counter++] = (*line);
			(*partition_config.node_in_current_block)[(*line)[0]-1] = 1; // line[0] == global_node_id
			delayed_lines_queue.pop_front();
		}
	}

	delete lines;
	return input;
}



#endif /*GRAPHIOSTREAM_H_*/
