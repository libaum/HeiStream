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

#include "definitions.h"
#include "priority_queue_interface.h"

const unsigned INVALID_POS = std::numeric_limits<unsigned>::max();
const float INVALID_SCORE = std::numeric_limits<float>::max();

class PQItem {
public:
    float buffer_score;
    int num_adj_partitioned;
     std::vector<LongNodeID> line;

    PQItem() : buffer_score(INVALID_SCORE) {}

    PQItem(float score, const std::vector<LongNodeID> &l, int adj_part)
        : buffer_score(score), num_adj_partitioned(adj_part), line(l) {
    }

    //  ~PQItem() {
    //      if (line != nullptr)
    //          delete line;
    //  }
    bool is_valid() const {
        return buffer_score != INVALID_SCORE;
    }

    std::vector<LongNodeID> &get_line() {
        return line;
    }

    void make_invalid() {
        buffer_score = INVALID_SCORE;
    }

    // void clean() { line.reset(); valid = false; }

    //  void clean() {
    //      if (line != nullptr) {
    //          delete line;
    //          line = nullptr;
    //      }
    //  }
};

class bucket_pq {
public:
    bucket_pq(const EdgeWeight &gain_span, LongNodeID max_node_id, LongNodeID max_pq_size);

    ~bucket_pq() = default;

    LongNodeID size();
    void insert(LongNodeID id, unsigned buffer_score);
    bool empty();

    unsigned maxValue();
    LongNodeID maxElement();
    LongNodeID deleteMax();

    void decreaseKey(LongNodeID node, unsigned newGain);
    void increaseKey(LongNodeID node, unsigned newGain);

    void changeKey(LongNodeID element, unsigned newKey);
    unsigned getKey(LongNodeID element);
    void deleteNode(LongNodeID node);

    bool contains(LongNodeID node);

private:
    LongNodeID m_elements;
    EdgeWeight m_gain_span;
    unsigned m_max_idx; // points to the non-empty bucket with the largest gain

    static constexpr LongNodeID VECTOR_THRESHOLD = 10000000; // Adjust based on performance tests
    bool use_vector;                                         // Decides whether to use vector or unordered_map

    // Hybrid data structure
    std::vector<std::pair<NodeID, unsigned>> m_queue_index_vec;
    std::unordered_map<LongNodeID, std::pair<NodeID, unsigned>> m_queue_index_map;
    std::vector<std::vector<LongNodeID>> m_buckets;
};

inline bucket_pq::bucket_pq(const EdgeWeight &buffer_score_span_input, LongNodeID num_nodes, LongNodeID max_pq_size)
    : m_elements(0), m_gain_span(buffer_score_span_input), m_max_idx(0) {

    use_vector = (num_nodes < VECTOR_THRESHOLD);

    if (use_vector) {
        m_queue_index_vec.resize(num_nodes, std::make_pair(INVALID_POS, 0));
    } else {
        m_queue_index_map.reserve(max_pq_size);
        m_queue_index_map.max_load_factor(0.7f);
    }
    m_buckets.resize(m_gain_span);
}

inline LongNodeID bucket_pq::size() {
    return m_elements;
}

inline void bucket_pq::insert(LongNodeID node, unsigned buffer_score) {
    unsigned address = buffer_score;
    assert(0 <= address && address < m_gain_span);

    if (address > m_max_idx) {
        m_max_idx = address;
    }

    m_buckets[address].push_back(node);
    if (use_vector) {
        m_queue_index_vec[node].first = m_buckets[address].size() - 1; // store position
        m_queue_index_vec[node].second = buffer_score;
    } else {
        m_queue_index_map.insert(
            std::make_pair(node, std::make_pair(m_buckets[address].size() - 1, buffer_score))
        );
    }

    m_elements++;
}

inline bool bucket_pq::empty() {
    return m_elements == 0;
}

inline unsigned bucket_pq::maxValue() {
    return m_max_idx;
}

inline LongNodeID bucket_pq::maxElement() {
    return m_buckets[m_max_idx].back();
}

inline LongNodeID bucket_pq::deleteMax() {
    LongNodeID node = m_buckets[m_max_idx].back();
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

inline void bucket_pq::decreaseKey(LongNodeID node, unsigned new_buffer_score) {
    changeKey(node, new_buffer_score);
}

inline void bucket_pq::increaseKey(LongNodeID node, unsigned new_buffer_score) {
    assert(0 <= new_buffer_score && new_buffer_score < m_gain_span);
    changeKey(node, new_buffer_score);
}

inline unsigned bucket_pq::getKey(LongNodeID node) {
    if (use_vector) {
        return m_queue_index_vec[node].second;
    } else {
        return m_queue_index_map[node].second;
    }
}

inline void bucket_pq::changeKey(LongNodeID node, unsigned new_buffer_score) {
    deleteNode(node);
    insert(node, new_buffer_score);
}

inline void bucket_pq::deleteNode(LongNodeID node) {
    if (use_vector) {
        LongNodeID vec_idx = node;
        assert(m_queue_index_vec[vec_idx].first != INVALID_POS);

        unsigned in_bucket_idx = m_queue_index_vec[vec_idx].first;
        unsigned old_buffer_score = m_queue_index_vec[vec_idx].second;
        unsigned address = old_buffer_score;

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

        unsigned in_bucket_idx = m_queue_index_map[node].first;
        unsigned old_buffer_score = m_queue_index_map[node].second;
        unsigned address = old_buffer_score;

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


inline bool bucket_pq::contains(LongNodeID node) {
    if (use_vector) {
        return m_queue_index_vec[node].first != INVALID_POS;
    } else {
        return m_queue_index_map.find(node) != m_queue_index_map.end();
    }
}

#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */