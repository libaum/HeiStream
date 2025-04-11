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
#include <sparsehash/dense_hash_map>


const float INVALID_SCORE = std::numeric_limits<float>::max();
const unsigned INVALID_POS = std::numeric_limits<unsigned>::max();
// const unsigned INVALID_POS_CACHED = std::numeric_limits<unsigned>::max() - 1;

inline unsigned discretize_score(float score) {
    // Use round instead of floor to handle precision better
    // return static_cast<unsigned>(std::round(score * config.bq_disc_factor));
    return static_cast<unsigned>(std::round(score * 100));
}

class PQItem {
public:
    float buffer_score;
    unsigned address;
    unsigned num_adj_partitioned;
    std::vector<LongNodeID> *adjacents;
    unsigned pos_in_bucket;

    // Default-Konstruktor
    PQItem() : buffer_score(INVALID_SCORE), address(0), num_adj_partitioned(0),
                adjacents(nullptr), pos_in_bucket(INVALID_POS) {}

    // PQItem(float score, const std::vector<LongNodeID> &l, unsigned adj_part, unsigned pos_in_bucket, unsigned address)
    //     : buffer_score(score), num_adj_partitioned(adj_part), adjacents(l), pos_in_bucket(pos_in_bucket), address(address) {
    //     adjacents = new std::vector<LongNodeID>(l);
    // }

    // Konstruktor für normale Referenzen (lvalue)
    PQItem(float score, std::vector<LongNodeID> &l, unsigned adj_part, unsigned pos_in_bucket, unsigned address)
        : buffer_score(score), address(address), num_adj_partitioned(adj_part), pos_in_bucket(pos_in_bucket) {
        // Übernehme Ownership des Vektors ohne Kopie
        // adjacents = new std::vector<LongNodeID>(std::move(l));
        adjacents = new std::vector<LongNodeID>(l);
    }

    // Move-Konstruktor für rvalue-Referenzen
    PQItem(float score, std::vector<LongNodeID>&& l, unsigned adj_part, unsigned pos_in_bucket, unsigned address)
        : buffer_score(score), address(address), num_adj_partitioned(adj_part), pos_in_bucket(pos_in_bucket) {
        // Übernehme Ownership des Vektors ohne Kopie
        adjacents = new std::vector<LongNodeID>(std::move(l));
    }

    // Copy-Konstruktor
    PQItem(const PQItem& other) : buffer_score(other.buffer_score), address(other.address), num_adj_partitioned(other.num_adj_partitioned), pos_in_bucket(other.pos_in_bucket) {
        // Kopiere Nachbarn, falls vorhanden
        adjacents = other.adjacents ? new std::vector<LongNodeID>(*other.adjacents) : nullptr;
    }

    // Move-Konstruktor
    PQItem(PQItem&& other) noexcept : buffer_score(other.buffer_score), address(other.address),
                num_adj_partitioned(other.num_adj_partitioned), adjacents(other.adjacents), pos_in_bucket(other.pos_in_bucket) {
        // Vermeidet doppeltes Freigeben
        other.adjacents = nullptr;
    }

    // Copy-Zuweisungsoperator
    PQItem& operator=(const PQItem& other) {
        if (this != &other) {
            buffer_score = other.buffer_score;
            address = other.address;
            num_adj_partitioned = other.num_adj_partitioned;
            pos_in_bucket = other.pos_in_bucket;

            // Alten Vector löschen
            delete adjacents;

            // Neuen Vector kopieren, falls vorhanden
            adjacents = other.adjacents ? new std::vector<LongNodeID>(*other.adjacents) : nullptr;
        }
        return *this;
    }

    // Move-Zuweisungsoperator
    PQItem& operator=(PQItem&& other) noexcept {
        if (this != &other) {
            buffer_score = other.buffer_score;
            address = other.address;
            num_adj_partitioned = other.num_adj_partitioned;
            pos_in_bucket = other.pos_in_bucket;

            // Alten Vector löschen
            delete adjacents;

            // Pointer übernehmen
            adjacents = other.adjacents;
            other.adjacents = nullptr;
        }
        return *this;
    }

    // Destruktor bleibt ähnlich, aber jetzt mit nullptr-Check
    ~PQItem() {
        if (adjacents != nullptr) {
            delete adjacents;  // Sicher für nullptr
        }
    }

    std::vector<LongNodeID>& get_adjacents() {
        return *adjacents;
    }

    unsigned get_address() const {
        return address;
    }

    // unsigned get_disc_buffer_score() const {
    //     return discretize_score(buffer_score);
    // }

    unsigned get_pos_in_bucket() const {
        return pos_in_bucket;
    }

    void set_pos_in_bucket(unsigned pos) {
        pos_in_bucket = pos;
    }

    void set_buffer_score(float new_buffer_score, unsigned new_address) {
        buffer_score = new_buffer_score;
        address = new_address;
    }

    // bool is_cached() const {
    //     return pos_in_bucket == INVALID_POS_CACHED;
    // }

    // bool is_valid() const {
    //     return pos_in_bucket != INVALID_POS_CACHED;
    // }

    // void set_cached() {
    //     pos_in_bucket = INVALID_POS_CACHED;
    // }


    // void clean() { adjacents.reset(); valid = false; }

    //  void clean() {
    //      if (adjacents != nullptr) {
    //          delete adjacents;
    //          adjacents = nullptr;
    //      }
    //  }
};

class bucket_pq {
public:
    bucket_pq(const EdgeWeight &gain_span, LongNodeID max_node_id, LongNodeID max_pq_size);

    ~bucket_pq() = default;

    LongNodeID size() const;
    void update_insert(LongNodeID node, float buffer_score);
    void insert(LongNodeID node, PQItem&& item);
    void insert(LongNodeID node, float buffer_score,  std::vector<LongNodeID>& adjacents, unsigned num_adj_partitioned);
    bool empty() const;

    unsigned maxValue() const;
    LongNodeID maxElement();
    LongNodeID deleteMax();

    void decreaseKey(LongNodeID node, unsigned newGain);
    void increaseKey(LongNodeID node, float newGain);

    void changeKey(LongNodeID element, float newKey);
    unsigned getKey(LongNodeID element);
    void deleteNode(LongNodeID node);

    bool contains(LongNodeID node);
    PQItem& getBufferItem(LongNodeID node);
    void completely_remove_node(LongNodeID node);

private:
    LongNodeID m_elements;
    EdgeWeight m_gain_span;
    unsigned m_max_idx; // points to the non-empty bucket with the largest gain

    std::unordered_map<LongNodeID, PQItem> m_queue_index_map;
    std::vector<std::vector<LongNodeID>> m_buckets;
};

inline bucket_pq::bucket_pq(const EdgeWeight &buffer_score_span_input, LongNodeID num_nodes, LongNodeID max_pq_size)
    : m_elements(0), m_gain_span(buffer_score_span_input), m_max_idx(0) {

    m_queue_index_map.reserve(max_pq_size);
    m_queue_index_map.max_load_factor(0.7f);
    m_buckets.resize(m_gain_span);
}

inline LongNodeID bucket_pq::size() const {
    return m_elements;
}


inline bool bucket_pq::contains(LongNodeID node) {
    return m_queue_index_map.find(node) != m_queue_index_map.end();
}


inline PQItem& bucket_pq::getBufferItem(LongNodeID node) {
    return m_queue_index_map[node];
}

inline void bucket_pq::completely_remove_node(LongNodeID node) {
    m_queue_index_map.erase(node);
}

inline void bucket_pq::update_insert(LongNodeID node, float buffer_score) {


    unsigned new_address = discretize_score(buffer_score);
    assert(0 <= new_address && new_address < m_gain_span);

    if (new_address > m_max_idx) {
        m_max_idx = new_address;
    }

    m_buckets[new_address].push_back(node);
    m_queue_index_map[node].set_buffer_score(buffer_score, new_address);
    m_queue_index_map[node].set_pos_in_bucket(m_buckets[new_address].size() - 1);

    m_elements++;
}


inline void bucket_pq::insert(LongNodeID node, float buffer_score,  std::vector<LongNodeID>& adjacents, unsigned num_adj_partitioned) {

    unsigned address = discretize_score(buffer_score);

    if (address > m_max_idx) {
        m_max_idx = address;
    }

    m_buckets[address].push_back(node);
    m_queue_index_map.emplace(
        node,
        PQItem(buffer_score, adjacents, num_adj_partitioned, m_buckets[address].size() - 1, address)
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

inline void bucket_pq::decreaseKey(LongNodeID node, unsigned new_buffer_score) {
    changeKey(node, new_buffer_score);
}

inline void bucket_pq::increaseKey(LongNodeID node, float new_buffer_score) {
    changeKey(node, new_buffer_score);
}

inline unsigned bucket_pq::getKey(LongNodeID node) {
    return m_queue_index_map[node].pos_in_bucket;
}

inline void bucket_pq::changeKey(LongNodeID node, float new_buffer_score) {
    deleteNode(node);
    update_insert(node, new_buffer_score);
}

inline void bucket_pq::deleteNode(LongNodeID node) {
    PQItem& buffer_item = m_queue_index_map[node];

    unsigned& in_bucket_idx = buffer_item.pos_in_bucket;
    unsigned& address = buffer_item.address;

    if (m_buckets[address].size() > 1) {
        // swap current element with last element and pop_back
        m_queue_index_map[m_buckets[address].back()].pos_in_bucket = in_bucket_idx; // update helper structure
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
    // buffer_item.set_cached();
}

#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */