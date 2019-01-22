#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "sv.h"

class Cluster {
    public:
    Cluster() = default;
    Cluster(const SV &_sv);
    SV cluster_represent;
    std::vector<SV> SVs;
    
    void merge(const SV &_sv);
    void merge(const Cluster &other);

    bool intersect(const Cluster &other, const float max_diff,
        int32_t max_dist) const;
};

#endif // CLUSTER_H