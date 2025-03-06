/******************************************************************************
 * bucket_pq.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

 #ifndef BUCKET_PQ_EM8YJPA9
 #define BUCKET_PQ_EM8YJPA9

 #include <limits>
 #include <unordered_map>
 #include <utility>

 #include "priority_queue_interface.h"

 const Count INVALID_POS = std::numeric_limits<Count>::max();
 const float INVALID_SCORE = std::numeric_limits<float>::max();

 class PQItem {
 public:
     float buffer_score;
     int num_adj_partitioned;
     std::vector<LongNodeID> *line;

     PQItem()
         : buffer_score(INVALID_SCORE),
           num_adj_partitioned(0),
           line(nullptr) {}

     ~PQItem() {
         if (line != nullptr)
             delete line;
     }

     void make_invalid() {
         buffer_score = INVALID_SCORE;
     }

     bool is_valid() const {
         return buffer_score != INVALID_SCORE;
         // return line != nullptr;
     }

     void clean() {
         if (line != nullptr) {
             delete line;
             line = nullptr;
         }
     }
 };

 class bucket_pq : public priority_queue_interface {
 public:
     bucket_pq(const EdgeWeight &gain_span, NodeID max_node_id);

     ~bucket_pq() override = default;

     NodeID size() override;
     void insert(NodeID id, Gain gain) override;
     bool empty() override;

     Gain maxValue() override;
     NodeID maxElement() override;
     NodeID deleteMax() override;

     void decreaseKey(NodeID node, Gain newGain) override;
     void increaseKey(NodeID node, Gain newGain) override;

     void changeKey(NodeID element, Gain newKey) override;
     Gain getKey(NodeID element) override;
     void deleteNode(NodeID node) override;

     bool contains(NodeID node) override;

 private:
     NodeID m_elements;
     EdgeWeight m_gain_span;
     unsigned m_max_idx; // points to the non-empty bucket with the largest gain

     static constexpr NodeID VECTOR_THRESHOLD = 10000000; // Adjust based on performance tests
     bool use_vector;                                    // Decides whether to use vector or unordered_map

     // Hybrid data structure
     std::vector<std::pair<Count, Gain>> m_queue_index_vec;
     std::unordered_map<NodeID, std::pair<Count, Gain>> m_queue_index_map;
     std::vector<std::vector<NodeID>> m_buckets;
 };

 inline bucket_pq::bucket_pq(const EdgeWeight &gain_span_input, NodeID num_nodes)
     : m_elements(0), m_gain_span(gain_span_input), m_max_idx(0) {


     use_vector = (num_nodes < VECTOR_THRESHOLD);

     if (use_vector) {
         m_queue_index_vec.resize(num_nodes, std::make_pair(INVALID_POS, 0));
     }
     m_buckets.resize(m_gain_span);
 }

 inline NodeID bucket_pq::size() {
     return m_elements;
 }

 inline void bucket_pq::insert(NodeID node, Gain gain) {
     unsigned address = gain;
     assert(0 <= address && address < m_gain_span);

     if (address > m_max_idx) {
         m_max_idx = address;
     }

     m_buckets[address].push_back(node);
     if (use_vector) {
         m_queue_index_vec[node].first = m_buckets[address].size() - 1; // store position
         m_queue_index_vec[node].second = gain;
     } else {
         m_queue_index_map[node].first = m_buckets[address].size() - 1; // store position
         m_queue_index_map[node].second = gain;
     }

     m_elements++;
 }

 inline bool bucket_pq::empty() {
     return m_elements == 0;
 }

 inline Gain bucket_pq::maxValue() {
     return m_max_idx;
 }

 inline NodeID bucket_pq::maxElement() {
     return m_buckets[m_max_idx].back();
 }

 inline NodeID bucket_pq::deleteMax() {
     NodeID node = m_buckets[m_max_idx].back();
     m_buckets[m_max_idx].pop_back();
     if (use_vector) {
         m_queue_index_vec[node].first = INVALID_POS;
     } else {
         m_queue_index_map.erase(node);
     }

     if (m_buckets[m_max_idx].size() == 0) {
         // update max_idx
         while (m_max_idx != 0) {
             m_max_idx--;
             if (m_buckets[m_max_idx].size() > 0) {
                 break;
             }
         }
     }

     m_elements--;
     return node;
 }

 inline void bucket_pq::decreaseKey(NodeID node, Gain new_gain) {
     changeKey(node, new_gain);
 }

 inline void bucket_pq::increaseKey(NodeID node, Gain new_gain) {
     assert(0 <= new_gain && new_gain <= m_gain_span);
     changeKey(node, new_gain);
 }

 inline Gain bucket_pq::getKey(NodeID node) {
     if (use_vector) {
         return m_queue_index_vec[node].second;
     } else {
         return m_queue_index_map[node].second;
     }
 }

 inline void bucket_pq::changeKey(NodeID node, Gain new_gain) {
     deleteNode(node);
     insert(node, new_gain);
 }

 inline void bucket_pq::deleteNode(NodeID node) {
     if (use_vector) {
         unsigned vec_idx = node;
         assert(m_queue_index_vec[vec_idx].first != INVALID_POS);

         Count in_bucket_idx = m_queue_index_vec[vec_idx].first;
         Gain old_gain = m_queue_index_vec[vec_idx].second;
         unsigned address = old_gain;

         if (m_buckets[address].size() > 1) {
             // swap current element with last element and pop_back
             m_queue_index_vec[m_buckets[address].back()].first = in_bucket_idx; // update helper structure
             std::swap(m_buckets[address][in_bucket_idx], m_buckets[address].back());
             m_buckets[address].pop_back();
         } else {
             // size is 1
             m_buckets[address].pop_back();
             if (address == m_max_idx) {
                 // update max_idx
                 while (m_max_idx != 0) {
                     m_max_idx--;
                     if (m_buckets[m_max_idx].size() > 0) {
                         break;
                     }
                 }
             }
         }

         m_elements--;
         m_queue_index_vec[vec_idx].first = INVALID_POS;

     } else {
         assert(m_queue_index_map.find(node) != m_queue_index_map.end());

         Count in_bucket_idx = m_queue_index_map[node].first;
         Gain old_gain = m_queue_index_map[node].second;
         unsigned address = old_gain;

         if (m_buckets[address].size() > 1) {
             // swap current element with last element and pop_back
             m_queue_index_map[m_buckets[address].back()].first = in_bucket_idx; // update helper structure
             std::swap(m_buckets[address][in_bucket_idx], m_buckets[address].back());
             m_buckets[address].pop_back();
         } else {
             // size is 1
             m_buckets[address].pop_back();
             if (address == m_max_idx) {
                 // update max_idx
                 while (m_max_idx != 0) {
                     m_max_idx--;
                     if (m_buckets[m_max_idx].size() > 0) {
                         break;
                     }
                 }
             }
         }

         m_elements--;
         m_queue_index_map.erase(node);
     }
 }

 inline bool bucket_pq::contains(NodeID node) {
     if (use_vector) {
         return m_queue_index_vec[node].first != INVALID_POS;
     } else {
         return m_queue_index_map.find(node) != m_queue_index_map.end();
     }
 }

 #endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */
