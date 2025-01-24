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
#include <list>
#include <cassert>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
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
#include "data_structure/priority_queues/delayed_lines_pq.h"

#define MIN(A,B) ((A)<(B))?(A):(B)
#define MAX(A,B) ((A)>(B))?(A):(B)

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

// CUTTANA HYPERPARAMETERS
const double THETA = 2;
const int D_MAX = 1000;
LongNodeID MAX_QUEUE_SIZE = 1000000;

const LongNodeID INVALID_NODE = std::numeric_limits<LongNodeID>::max();

// Constants for bucket queue
const float MAX_BUFFER_SCORE = 3.0f;
const int DISCRETIZATION_FACTOR = 100; // Two decimal places
const int MAX_BUCKET_INDEX = static_cast<int>(std::floor(MAX_BUFFER_SCORE * DISCRETIZATION_FACTOR));

// BucketQueue class definition
class BucketQueue {
public:
    BucketQueue()
        : buckets(MAX_BUCKET_INDEX + 1), current_bucket(MAX_BUCKET_INDEX), size_(0) {}

	void push(LongNodeID node, float score) {
		int bucket_idx = discretize_score(score);

		// Clamp the bucket index to valid range
		bucket_idx = std::max(0, std::min(bucket_idx, max_bucket_index));

		// Push the node to the corresponding bucket
		buckets[bucket_idx].push_back(node);

		// Update current_bucket if necessary
		current_bucket = std::max(current_bucket, bucket_idx);
		size_++;
	}

    bool empty() const {
        return size_ == 0;
    }

    LongNodeID pop() {
        while (current_bucket >= 0) {
            auto& current_bucket_list = buckets[current_bucket];

            if (!current_bucket_list.empty()) {
                LongNodeID node = current_bucket_list.back();
                current_bucket_list.pop_back();
                size_--;

                // Decrease current_bucket if the current bucket is empty
                if (buckets[current_bucket].empty()) {
                    --current_bucket;
                }

                return node;
            }

            --current_bucket;
        }

        return INVALID_NODE; // Indicates the queue is empty
    }

    size_t size() const {
        return size_;
    }

private:
    std::vector<std::list<LongNodeID>> buckets;
    int current_bucket;
    int max_bucket_index = MAX_BUCKET_INDEX;
    size_t size_;

    int discretize_score(float score) const {
        // Use round instead of floor to handle precision better
        return static_cast<int>(std::round(score * DISCRETIZATION_FACTOR));
    }
};


void config_multibfs_initial_partitioning(PartitionConfig & partition_config);


float calc_buffer_score(PartitionConfig &partition_config, std::vector<LongNodeID> &cur_line) {
	int global_node_id = cur_line[0];
	float degree = cur_line.size() - 1;
	float cnt_adj_partitioned = 0;
	bool adj_is_partitioned;
	for (LongNodeID &global_adj_id : cur_line) {
		if (global_adj_id == global_node_id) continue;
		adj_is_partitioned = (*partition_config.stream_nodes_assign)[global_adj_id-1] !=  INVALID_PARTITION;
		if(adj_is_partitioned) {
			cnt_adj_partitioned++;
		}
	}
	float buffer_score = degree / D_MAX + THETA * cnt_adj_partitioned / degree; // Range: [0, 3] (first term: [0, 1], second term: [0, THETA])
	return buffer_score;
}

// Update the priority value of the neighbours of the node that was just partitioned in the priority queue
void update_neighbours_priority(PartitionConfig &partition_config, std::vector<LongNodeID> &line, std::vector<float> &node_id_to_buffer_score, BucketQueue &pq, std::vector<std::vector<LongNodeID>> &node_id_to_line) {
	LongNodeID global_node_id = line[0];
	int degree = line.size() - 1;
	for (auto it = line.begin() + 1; it != line.end(); ++it) {
		LongNodeID adj_id = *it;

		// Update buffer score of neighbours
		bool is_partitioned = (*partition_config.stream_nodes_assign)[adj_id-1] != INVALID_PARTITION;
		bool is_in_pq = node_id_to_buffer_score[adj_id] != UNDEFINED_LONGNODE;
		if (!is_partitioned && is_in_pq) {
			int adj_degree = node_id_to_line[adj_id].size() - 1;
			float updated_buffer_score = node_id_to_buffer_score[adj_id] + THETA / adj_degree ;
			// if (std::fabs(updated_buffer_score - node_id_to_buffer_score[adj_id]) > 1e-4) { // Threshold to avoid minor changes
			// 	node_id_to_buffer_score[adj_id] = updated_buffer_score;
			// 	pq.push(adj_id, updated_buffer_score);
			// }
			node_id_to_buffer_score[adj_id] = updated_buffer_score;
			pq.push(adj_id, updated_buffer_score);

		}

	}
}

void partition_node(PartitionConfig &partition_config, std::vector<LongNodeID> &line) {
	LongNodeID global_node_id = line[0];
	std::vector<int> hash_map(partition_config.k, 0);

	for (auto it = line.begin() + 1; it != line.end(); ++it) {
		LongNodeID adj_id = *it;
		PartitionID adj_part = (*partition_config.stream_nodes_assign)[adj_id-1];
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
		(*partition_config.stream_nodes_assign)[global_node_id-1] = best_partition;

		// Step 5: Update partition load
		(*partition_config.stream_blocks_weight)[best_partition]++;

	}

}



int main(int argn, char **argv) {
        PartitionConfig partition_config;
        std::string graph_filename;
        timer t, processing_t, io_t;
	EdgeWeight total_edge_cut = 0;
        double global_mapping_time = 0;
	double buffer_mapping_time = 0;
	double buffer_io_time = 0;
        quality_metrics qm;
	balance_configuration bc;
	// std::vector<std::vector<LongNodeID>>* input = nullptr;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;

        int ret_code = parse_parameters(argn, argv,
                                        partition_config,
                                        graph_filename,
                                        is_graph_weighted,
                                        suppress_output, recursive);

        if(ret_code) {
                return 0;
        }

        std::ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
                std::cout.rdbuf(ofs.rdbuf());
        }
	srand(partition_config.seed);
	random_functions::setSeed(partition_config.seed);

        partition_config.LogDump(stdout);
	partition_config.stream_input = true;
	graph_access *G = new graph_access();

	timer first_pass_t, second_pass_t, updating_adj_t, partitioning_t, calc_buffer_score_t;
	double first_pass_time = 0;
	double second_pass_time = 0;
	double updating_adj_time = 0;
	double partitioning_time = 0;
	double calc_buffer_score_time = 0;

	int &passes = partition_config.num_streams_passes;
	for (partition_config.restream_number=0; partition_config.restream_number<passes; partition_config.restream_number++) {

		// ***************************** IO operations ***************************************
		io_t.restart();
		graph_io_stream::readFirstLineStream(partition_config, graph_filename, total_edge_cut);


		double avg_block_size = static_cast<double>(partition_config.number_of_nodes) / partition_config.k;
		partition_config.max_block_weight = static_cast<int>(std::ceil((1.0 + partition_config.imbalance / 100) * avg_block_size));

		buffer_io_time += io_t.elapsed();
    	// std::priority_queue<lines_pq_node, std::vector<lines_pq_node>, CompareBufferScore> pq_delayed_lines;
    	BucketQueue bucket_pq;
    	std::vector<float> node_id_to_buffer_score(partition_config.number_of_nodes, UNDEFINED_LONGNODE); // if -INFINITY: not in priority queue, else: index in priority queue with corresponding priority
		std::vector<std::vector<LongNodeID>> node_id_to_line(partition_config.number_of_nodes, std::vector<LongNodeID>(0));
		// LongNodeID queue_size = 0;

		LongNodeID node_counter = 0;
		std::unique_ptr<buffered_input> ss2 = nullptr;
		std::vector<LongNodeID> cur_line;

		auto lines = std::make_unique<std::vector<std::string>>(1);

		int processed_nodes = 0;
		first_pass_t.restart();
		while (partition_config.remaining_stream_nodes) {
			// Load a line from the stream
			std::getline(*(partition_config.stream_in),(*lines)[0]);
			if ((*lines)[0][0] == '%') { // a comment in the file
				continue;
			}
			partition_config.remaining_stream_nodes--;
			LongNodeID global_node_id = ++partition_config.total_nodes_loaded;

			ss2 = std::make_unique<buffered_input>(lines.get());
			ss2->simple_scan_line(cur_line);
			cur_line.insert(cur_line.begin(), global_node_id);

			int degree = cur_line.size() - 1;
			if (degree > D_MAX) {
				// std::cout << "degree > D_MAX. Partitioning node directly, degree: " << degree << std::endl;
				// Partition node directly
				partitioning_t.restart();
				partition_node(partition_config, cur_line);
				partitioning_time += partitioning_t.elapsed();
				updating_adj_t.restart();
				update_neighbours_priority(partition_config, cur_line, node_id_to_buffer_score, bucket_pq, node_id_to_line);
				updating_adj_time += updating_adj_t.elapsed();

				node_id_to_buffer_score[global_node_id] = UNDEFINED_LONGNODE;
				node_id_to_line[global_node_id].clear();
				continue;
			} else if (bucket_pq.size() >= MAX_QUEUE_SIZE) {
				// Make space for new node
				// Remove node from queue by popping
				LongNodeID node_id_to_remove = bucket_pq.pop();

                // Check if the node is already partitioned
				bool is_not_partitioned_yet = (*partition_config.stream_nodes_assign)[node_id_to_remove - 1] == INVALID_PARTITION;
				if (is_not_partitioned_yet) {
					// Partition the node
					partitioning_t.restart();
					partition_node(partition_config, node_id_to_line[node_id_to_remove]);
					partitioning_time += partitioning_t.elapsed();

					updating_adj_t.restart();
					update_neighbours_priority(partition_config, node_id_to_line[node_id_to_remove], node_id_to_buffer_score, bucket_pq, node_id_to_line);
					updating_adj_time += updating_adj_t.elapsed();

					node_id_to_buffer_score[node_id_to_remove] = UNDEFINED_LONGNODE;
					node_id_to_line[node_id_to_remove].clear();
					// queue_size--;
				}


			}

			// Calculate priority of the node and push into BucketQueue
			calc_buffer_score_t.restart();
			float buffer_score = calc_buffer_score(partition_config, cur_line);
			calc_buffer_score_time += calc_buffer_score_t.elapsed();
			// std::cout << "Buffer score: " << buffer_score << std::endl;
			// pq_delayed_lines.push(lines_pq_node(cur_line, buffer_score));
			bucket_pq.push(global_node_id, buffer_score);
			node_id_to_buffer_score[global_node_id] = buffer_score;
			node_id_to_line[global_node_id] = cur_line;
			// queue_size++;

            // Optionally print progress
			processed_nodes++;
			if (processed_nodes % 100000 == 0) {
				// std::cout << "Processed " << processed_nodes << " nodes." << std::endl;
			}
		}
		cur_line.clear();
		first_pass_time += first_pass_t.elapsed();

		// std::cout << "Finished loading all nodes. Emptying queue now. Size: " << bucket_pq.size() << std::endl;// << ", queue size: " << queue_size << std::endl;
		second_pass_t.restart();
		while (!bucket_pq.empty()) {// && queue_size != 0) {
			LongNodeID node_id_to_remove = bucket_pq.pop();

			bool is_already_partitioned = (*partition_config.stream_nodes_assign)[node_id_to_remove - 1] != INVALID_PARTITION;
            if (is_already_partitioned) {
                continue;
            }

			// Partition the node
            partitioning_t.restart();
            partition_node(partition_config, node_id_to_line[node_id_to_remove]);
            partitioning_time += partitioning_t.elapsed();

            // Update neighbors' buffer scores
            updating_adj_t.restart();
            update_neighbours_priority(partition_config, node_id_to_line[node_id_to_remove], node_id_to_buffer_score, bucket_pq, node_id_to_line);
            updating_adj_time += updating_adj_t.elapsed();

            // queue_size--;

			// Optionally print queue size
			if (bucket_pq.size() % 100000 == 0) {
				// std::cout << "Queue size: " << bucket_pq.size() << std::endl;
			}
		}
		second_pass_time += second_pass_t.elapsed();

	// 	std::cout << "Finished emptying queue. Size: " << bucket_pq.size() << std::endl;
	}
	// std::cout << "First pass time: " << first_pass_time << std::endl;
	// std::cout << "Second pass time: " << second_pass_time << std::endl;
	// std::cout << "Updating adj time: " << updating_adj_time << std::endl;
	// std::cout << "Partitioning time: " << partitioning_time << std::endl;
	// std::cout << "Calc buffer score time: " << calc_buffer_score_time << std::endl;
	double total_time = processing_t.elapsed();
	long maxRSS = getMaxRSS();

	// Check if all nodes are assigned
	for (LongNodeID i = 0; i < partition_config.number_of_nodes; i++) {
		if ((*partition_config.stream_nodes_assign)[i] == INVALID_PARTITION) {
			std::cout << "Node " << i << " is not assigned." << std::endl;
		}
	}

	graph_io_stream::streamEvaluatePartition(partition_config, graph_filename, total_edge_cut);

		std::cout << total_time << " " << total_edge_cut;
		if (maxRSS != -1) {
			std::cout << " " << maxRSS << " " << qm.balance_full_stream(*partition_config.stream_blocks_weight) << std::endl;
		} else {
			std::cout << std::endl;
		}
        // write the partition to the disc
        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmppartition" << partition_config.k;
        } else {
                filename << partition_config.filename_output;
        }

        if (!partition_config.suppress_output) {
                // graph_io_stream::writePartitionStream(partition_config, filename.str());
        } else {
                std::cout << "No partition will be written as output." << std::endl;
        }

	delete G;
	// delete bucket_pq;

	if (partition_config.add_blocks_weight != nullptr) {
		delete partition_config.add_blocks_weight;
	}
	if (partition_config.node_in_current_block != nullptr) {
		delete partition_config.node_in_current_block;
	}
	if (partition_config.stream_nodes_assign != nullptr) {
		delete partition_config.stream_nodes_assign;
	}
	if (partition_config.local_to_global_map != nullptr) {
		delete partition_config.local_to_global_map;
	}
	if (partition_config.ghostkey_to_edges != nullptr) {
		delete partition_config.ghostkey_to_edges;
	}
	if (partition_config.stream_blocks_weight != nullptr) {
		delete partition_config.stream_blocks_weight;
	}
	if (partition_config.stream_in != nullptr) {
		delete partition_config.stream_in;
	}



	return 0;
}

void config_multibfs_initial_partitioning(PartitionConfig & partition_config) {
	if (partition_config.initial_part_multi_bfs && partition_config.curr_batch >= 2) {
		partition_config.initial_partitioning_type = INITIAL_PARTITIONING_MULTIBFS;
	}
}


