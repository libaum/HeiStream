#ifndef DELAYED_LINES_PQ
#define DELAYED_LINES_PQ

#include "lib/definitions.h"

struct lines_pq_node {
    std::vector<LongNodeID> line;
    double buffer_score;

    lines_pq_node(std::vector<LongNodeID> line, double buffer_score) : line(line), buffer_score(buffer_score) {}
};

struct CompareBufferScore {
    bool operator()(lines_pq_node const &n1, lines_pq_node const &n2) {
        // Higher priority first
        return n1.buffer_score < n2.buffer_score;
    }
};

#endif /* end of include guard: DELAYED_LINES_PQ */