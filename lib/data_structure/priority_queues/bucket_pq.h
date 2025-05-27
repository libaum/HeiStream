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

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>

#include "definitions.h"
#include "priority_queue_interface.h"
#include <sparsehash/dense_hash_map>


// #define ENABLE_TIMING  // Kommentiere diese Zeile aus, um Timing zu deaktivieren

#ifdef ENABLE_TIMING
#define START_TIMER(stats_ref) ScopedTimer timer(stats_ref)
#else
#define START_TIMER(stats_ref)
#endif


const float INVALID_SCORE = std::numeric_limits<float>::max();
const unsigned INVALID_POS = std::numeric_limits<unsigned>::max();
// const unsigned INVALID_POS_CACHED = std::numeric_limits<unsigned>::max() - 1;


class PQItem {
public:
    float buffer_score;
    unsigned bucket_idx;
    unsigned num_adj_partitioned;
    std::vector<LongNodeID> *adjacents;
    unsigned pos_in_bucket;

    // Default-Konstruktor
    PQItem() : buffer_score(INVALID_SCORE), bucket_idx(0), num_adj_partitioned(0),
                adjacents(nullptr), pos_in_bucket(INVALID_POS) {}

    // Konstruktor für normale Referenzen (lvalue)
    PQItem(float score, std::vector<LongNodeID> &l, unsigned adj_part, unsigned pos_in_bucket, unsigned bucket_idx)
        : buffer_score(score), bucket_idx(bucket_idx), num_adj_partitioned(adj_part), pos_in_bucket(pos_in_bucket) {
        adjacents = new std::vector<LongNodeID>(l);
    }

    // Move-Konstruktor für rvalue-Referenzen
    PQItem(float score, std::vector<LongNodeID>&& l, unsigned adj_part, unsigned pos_in_bucket, unsigned bucket_idx)
        : buffer_score(score), bucket_idx(bucket_idx), num_adj_partitioned(adj_part), pos_in_bucket(pos_in_bucket) {
        // Übernehme Ownership des Vektors ohne Kopie
        adjacents = new std::vector<LongNodeID>(std::move(l));
    }

    // Copy-Konstruktor
    PQItem(const PQItem& other) : buffer_score(other.buffer_score), bucket_idx(other.bucket_idx), num_adj_partitioned(other.num_adj_partitioned), pos_in_bucket(other.pos_in_bucket) {
        // Kopiere Nachbarn, falls vorhanden
        adjacents = other.adjacents ? new std::vector<LongNodeID>(*other.adjacents) : nullptr;
    }

    // Move-Konstruktor
    PQItem(PQItem&& other) noexcept : buffer_score(other.buffer_score), bucket_idx(other.bucket_idx),
                num_adj_partitioned(other.num_adj_partitioned), adjacents(other.adjacents), pos_in_bucket(other.pos_in_bucket) {
        // Vermeidet doppeltes Freigeben
        other.adjacents = nullptr;
    }

    // Copy-Zuweisungsoperator
    PQItem& operator=(const PQItem& other) {
        if (this != &other) {
            buffer_score = other.buffer_score;
            bucket_idx = other.bucket_idx;
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
            bucket_idx = other.bucket_idx;
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
        return bucket_idx;
    }

    unsigned get_pos_in_bucket() const {
        return pos_in_bucket;
    }

    void set_pos_in_bucket(unsigned pos) {
        pos_in_bucket = pos;
    }

    void set_buffer_score(float new_buffer_score, unsigned new_bucket_idx) {
        buffer_score = new_buffer_score;
        bucket_idx = new_bucket_idx;
    }
};

class bucket_pq {
public:
    bucket_pq(const EdgeWeight &gain_span, LongNodeID max_node_id, LongNodeID max_pq_size, unsigned bq_disc_factor);

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
    // bool contains_with_it(LongNodeID node, std::unordered_map<LongNodeID, PQItem>::iterator& it);
    PQItem& getBufferItem(LongNodeID node);
    void completely_remove_node(LongNodeID node);

    unsigned discretize_score(float score) const;

    void print_statistics() const;
    void reset_statistics();

    std::unordered_map<LongNodeID, PQItem>& get_queue_index_map() {
        return m_queue_index_map;
    }


private:
    LongNodeID m_elements;
    EdgeWeight m_gain_span;
    unsigned m_max_idx; // points to the non-empty bucket with the largest gain
    unsigned disc_factor;


    std::unordered_map<LongNodeID, PQItem> m_queue_index_map;
    std::vector<std::vector<LongNodeID>> m_buckets;

        // Timing-Strukturen
    struct OperationStats {
        double total_time = 0.0;       // Gesamtzeit in Sekunden
        uint64_t call_count = 0;       // Anzahl der Aufrufe
        double min_time = std::numeric_limits<double>::max();  // Minimale Zeit
        double max_time = 0.0;         // Maximale Zeit

        void update(double time_seconds) {
            total_time += time_seconds;
            call_count++;
            min_time = std::min(min_time, time_seconds);
            max_time = std::max(max_time, time_seconds);
        }

        double avg_time() const {
            return call_count > 0 ? total_time / call_count : 0.0;
        }
    };

    // Timer für verschiedene Operationen
    mutable struct {
        OperationStats insert;
        OperationStats update_insert;
        OperationStats delete_max;
        OperationStats delete_node;
        OperationStats change_key;
        OperationStats contains;
        OperationStats get_buffer_item;
    } stats;

    // Timer-Hilfsfunktion
    class ScopedTimer {
    private:
        OperationStats& stats;
        std::chrono::high_resolution_clock::time_point start;
    public:
        ScopedTimer(OperationStats& stats_ref) : stats(stats_ref),
                                            start(std::chrono::high_resolution_clock::now()) {}
        ~ScopedTimer() {
            #ifdef ENABLE_TIMING
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            stats.update(elapsed.count());
            #endif
        }
    };
};

inline bucket_pq::bucket_pq(const EdgeWeight &buffer_score_span_input, LongNodeID num_nodes, LongNodeID max_pq_size, unsigned bq_disc_factor)
    : m_elements(0), m_gain_span(buffer_score_span_input), m_max_idx(0), disc_factor(bq_disc_factor) {

    m_queue_index_map.reserve(max_pq_size);
    m_queue_index_map.max_load_factor(0.7f);
    m_buckets.resize(m_gain_span);

}

inline unsigned bucket_pq::discretize_score(float score) const{
    // Use round instead of floor to handle precision better
    // return static_cast<unsigned>(std::round(score * config.bq_disc_factor));
    return static_cast<unsigned>(std::round(score * disc_factor));
}

inline LongNodeID bucket_pq::size() const {
    return m_elements;
}


inline bool bucket_pq::contains(LongNodeID node) {
    START_TIMER(stats.contains);

    return m_queue_index_map.find(node) != m_queue_index_map.end();
}

// inline bool bucket_pq::contains_with_it(LongNodeID node, std::unordered_map<LongNodeID, PQItem>::iterator& it) {
//     START_TIMER(stats.contains);

//     it = m_queue_index_map.find(node);

//     return it != m_queue_index_map.end();
// }

inline PQItem& bucket_pq::getBufferItem(LongNodeID node) {
    START_TIMER(stats.get_buffer_item);

    return m_queue_index_map[node];
}

inline void bucket_pq::completely_remove_node(LongNodeID node) {
    m_queue_index_map.erase(node);
}

inline void bucket_pq::update_insert(LongNodeID node, float buffer_score) {

    START_TIMER(stats.update_insert);

    unsigned new_bucket_idx = discretize_score(buffer_score);
    assert(0 <= new_bucket_idx && new_bucket_idx < m_gain_span);

    if (new_bucket_idx > m_max_idx) {
        m_max_idx = new_bucket_idx;
    }

    m_buckets[new_bucket_idx].push_back(node);
    m_queue_index_map[node].set_buffer_score(buffer_score, new_bucket_idx);
    m_queue_index_map[node].set_pos_in_bucket(m_buckets[new_bucket_idx].size() - 1);

    m_elements++;
}


inline void bucket_pq::insert(LongNodeID node, float buffer_score,  std::vector<LongNodeID>& adjacents, unsigned num_adj_partitioned, unsigned bucket_idx) {
    START_TIMER(stats.insert);

    if (bucket_idx > m_max_idx) {
        m_max_idx = bucket_idx;
    }

    m_buckets[bucket_idx].push_back(node);
    m_queue_index_map.emplace(
        node,
        PQItem(buffer_score, adjacents, num_adj_partitioned, m_buckets[bucket_idx].size() - 1, bucket_idx)
    );

    // m_queue_index_map[node] = PQItem(buffer_score, adjacents, num_adj_partitioned, m_buckets[bucket_idx].size() - 1, bucket_idx);

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
    START_TIMER(stats.delete_max);

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
    return m_queue_index_map[node].pos_in_bucket;
}

inline void bucket_pq::changeKey(LongNodeID node, float new_buffer_score) {
    START_TIMER(stats.change_key);

    deleteNode(node);
    update_insert(node, new_buffer_score);
}

inline void bucket_pq::deleteNode(LongNodeID node) {
    START_TIMER(stats.delete_node);

    PQItem& buffer_item = m_queue_index_map[node];

    unsigned& in_bucket_idx = buffer_item.pos_in_bucket;
    unsigned& bucket_idx = buffer_item.bucket_idx;

    if (m_buckets[bucket_idx].size() > 1) {
        // Swap current element with last element and pop_back
        m_queue_index_map[m_buckets[bucket_idx].back()].pos_in_bucket = in_bucket_idx; // update helper structure
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




inline void bucket_pq::print_statistics() const {
#ifdef ENABLE_TIMING
    std::cout << "===== bucket_pq Laufzeitstatistiken =====" << std::endl;
    std::cout << "Aktueller PQ-Größe: " << m_elements << " Elemente" << std::endl;

    std::cout << "\nOperation          Aufrufe     Gesamt (s)    Durchschnitt (µs)    Min (µs)     Max (µs)" << std::endl;
    std::cout << "-----------------------------------------------------------------------------------------" << std::endl;

    auto print_stats = [](const char* name, const OperationStats& s) {
        std::cout << std::left << std::setw(18) << name
                  << std::right << std::setw(12) << s.call_count
                  << std::setw(14) << std::fixed << std::setprecision(6) << s.total_time
                  << std::setw(20) << std::fixed << std::setprecision(3) << s.avg_time() * 1e6
                  << std::setw(12) << std::fixed << std::setprecision(3) << s.min_time * 1e6
                  << std::setw(12) << std::fixed << std::setprecision(3) << s.max_time * 1e6
                  << std::endl;
    };

    print_stats("insert", stats.insert);
    print_stats("update_insert", stats.update_insert);
    print_stats("delete_max", stats.delete_max);
    print_stats("delete_node", stats.delete_node);
    print_stats("change_key", stats.change_key);
    print_stats("contains", stats.contains);
    print_stats("get_buffer_item", stats.get_buffer_item);

    // Gesamtstatistik
    double total_time = stats.insert.total_time + stats.update_insert.total_time +
                        stats.delete_max.total_time + stats.delete_node.total_time +
                        stats.change_key.total_time + stats.contains.total_time +
                        stats.get_buffer_item.total_time;

    uint64_t total_calls = stats.insert.call_count + stats.update_insert.call_count +
                         stats.delete_max.call_count + stats.delete_node.call_count +
                         stats.change_key.call_count + stats.contains.call_count +
                         stats.get_buffer_item.call_count;

    std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(18) << "GESAMT"
              << std::right << std::setw(12) << total_calls
              << std::setw(14) << std::fixed << std::setprecision(6) << total_time
              << std::endl;

    std::cout << "\nDurchschnittliche Zeit pro Element:" << std::endl;
    std::cout << " - Einfügen: " << (stats.insert.call_count > 0 ? stats.insert.total_time / stats.insert.call_count * 1e6 : 0.0) << " µs" << std::endl;
    std::cout << " - Entfernen: " << (stats.delete_max.call_count > 0 ? stats.delete_max.total_time / stats.delete_max.call_count * 1e6 : 0.0) << " µs" << std::endl;
#else
    std::cout << "Timer deactivated. Define ENABLE_TIMING for bucket PQ runtime stats." << std::endl;
#endif
}

inline void bucket_pq::reset_statistics() {
    stats = {}; // Zurücksetzen aller Statistiken auf Initialwerte
}

#endif /* end of include guard: BUCKET_PQ_EM8YJPA9 */