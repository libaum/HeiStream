/******************************************************************************
 * bucket_pq.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BUCKET_PQ_EM8YJPA9
#define BUCKET_PQ_EM8YJPA9

#include <memory>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>

#include "partition/partition_config.h"
#include "definitions.h"
#include "priority_queue_interface.h"
#include "pq_item.h"
#include "pq_storage_interface.h"
#include "vector_storage.h"
#include "map_storage.h"

class bucket_pq {
public:
    bucket_pq(PartitionConfig& partition_config, const unsigned &gain_span, LongNodeID max_node_id, LongNodeID max_buffer_size, unsigned bq_disc_factor);

    ~bucket_pq() = default;

    LongNodeID size() const;
    void update_insert(LongNodeID node, float buffer_score);
    void insert(LongNodeID node, float buffer_score,  std::vector<LongNodeID>& adjacents, unsigned num_adj_partitioned, unsigned bucket_idx);
    bool empty() const;

    unsigned maxValue() const;
    LongNodeID maxElement();
    LongNodeID deleteMax();

    void increaseKey(LongNodeID node, float newGain);

    void changeKey(LongNodeID element, float newKey);
    unsigned getKey(LongNodeID element);
    void deleteNode(LongNodeID node);

    bool contains(LongNodeID node);
    PQItem& getBufferItem(LongNodeID node);
    void completely_remove_node(LongNodeID node);

    unsigned discretize_score(float score) const;

    void print_statistics() const;
    void reset_statistics();

    std::unique_ptr<PQStorageBase<LongNodeID, PQItem>>& get_storage() {
        return m_storage;
    }

private:
    LongNodeID m_elements;
    unsigned m_gain_span;
    unsigned m_max_idx; // points to the non-empty bucket with the largest gain
    unsigned disc_factor;
    PartitionConfig &config;

    std::unique_ptr<PQStorageBase<LongNodeID, PQItem>> m_storage;
    std::vector<std::vector<LongNodeID>> m_buckets;
};

inline bucket_pq::bucket_pq(PartitionConfig& partition_config, const unsigned &buffer_score_span_input, LongNodeID num_nodes, LongNodeID max_buffer_size, unsigned bq_disc_factor)
    : config(partition_config), m_elements(0), m_gain_span(buffer_score_span_input), m_max_idx(0), disc_factor(bq_disc_factor) {

    // Storage-Typ zur Runtime wählen
    switch (config.bpq_storage_type) {
        // case BPQStorageType::BPQ_STORAGE_VECTOR:
        //     m_storage = std::make_unique<VectorStorage<LongNodeID, PQItem>>(num_nodes);
        //     break;
        case BPQStorageType::BPQ_STORAGE_UNORDERED_MAP:
        default:
            m_storage = std::make_unique<MapStorage<LongNodeID, PQItem>>();
            break;
    }

    m_storage->set_max_load_factor(0.7f);
    m_storage->reserve(max_buffer_size);
    m_buckets.resize(m_gain_span);

}

inline unsigned bucket_pq::discretize_score(float score) const{
    // Use round instead of floor to handle precision better
    return static_cast<unsigned>(std::round(score * disc_factor));
}

inline LongNodeID bucket_pq::size() const {
    return m_elements;
}


inline bool bucket_pq::contains(LongNodeID node) {
    return m_storage->contains(node);
}

inline PQItem& bucket_pq::getBufferItem(LongNodeID node) {
    return (*m_storage)[node];
}

inline void bucket_pq::completely_remove_node(LongNodeID node) {
    m_storage->erase(node);
}

inline void bucket_pq::update_insert(LongNodeID node, float buffer_score) {

    unsigned new_bucket_idx = std::min(discretize_score(buffer_score), m_gain_span - 1);
    assert(0 <= new_bucket_idx && new_bucket_idx < m_gain_span);

    if (new_bucket_idx > m_max_idx) {
        m_max_idx = new_bucket_idx;
    }

    m_buckets[new_bucket_idx].push_back(node);
    (*m_storage)[node].set_buffer_score(buffer_score, new_bucket_idx);
    (*m_storage)[node].set_pos_in_bucket(m_buckets[new_bucket_idx].size() - 1);

    m_elements++;
}


inline void bucket_pq::insert(LongNodeID node, float buffer_score,  std::vector<LongNodeID>& adjacents, unsigned num_adj_partitioned, unsigned bucket_idx) {
    if (bucket_idx > m_max_idx) {
        m_max_idx = bucket_idx;
    }

    m_buckets[bucket_idx].push_back(node);
    m_storage->emplace(
        node,
        PQItem(buffer_score, adjacents, num_adj_partitioned, m_buckets[bucket_idx].size() - 1, bucket_idx)
    );

    m_elements++;
}

inline bool bucket_pq::empty() const {
    return m_elements == 0;
}

inline unsigned bucket_pq::maxValue() const {
    return m_max_idx;
}

inline LongNodeID bucket_pq::maxElement() {
    return m_buckets[m_max_idx].back();
}

inline LongNodeID bucket_pq::deleteMax() {
    LongNodeID node = m_buckets[m_max_idx].back();
    m_buckets[m_max_idx].pop_back();

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

inline void bucket_pq::increaseKey(LongNodeID node, float new_buffer_score) {
    changeKey(node, new_buffer_score);
}

inline unsigned bucket_pq::getKey(LongNodeID node) {
    return (*m_storage)[node].pos_in_bucket;
}

inline void bucket_pq::changeKey(LongNodeID node, float new_buffer_score) {
    deleteNode(node);
    update_insert(node, new_buffer_score);
}

inline void bucket_pq::deleteNode(LongNodeID node) {
    PQItem& buffer_item = (*m_storage)[node];

    unsigned& in_bucket_idx = buffer_item.pos_in_bucket;
    unsigned& bucket_idx = buffer_item.bucket_idx;

    if (m_buckets[bucket_idx].size() > 1) {
        // Swap current element with last element and pop_back
        (*m_storage)[m_buckets[bucket_idx].back()].pos_in_bucket = in_bucket_idx; // update helper structure
        std::swap(m_buckets[bucket_idx][in_bucket_idx], m_buckets[bucket_idx].back());
        m_buckets[bucket_idx].pop_back();
    } else {
        // size is 1
        m_buckets[bucket_idx].pop_back();
        if (bucket_idx == m_max_idx) {
            // Update max_idx
            while (m_max_idx != 0) {
                m_max_idx--;
                if (m_buckets[m_max_idx].size() > 0) {
                    break;
                }
            }
        }
    }

    m_elements--;
}

#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */