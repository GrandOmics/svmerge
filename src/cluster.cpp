#include "cluster.h"
#include <iostream>

Cluster::Cluster(const SV &_sv) {
    cluster_represent = _sv;
    SVs.push_back(_sv);
}

bool Cluster::intersect(const Cluster &other,
    const float &max_diff, int32_t max_dist, const float &min_overlap) const
{
    if (SV_intersect(cluster_represent, other.cluster_represent,
        max_diff, max_dist, min_overlap))
    {
        return true;
    }
    return false;
}

void Cluster::merge(const SV &_sv) {
    SVs.push_back(_sv);
}

void Cluster::merge(const Cluster &other) {
    SVs.push_back(other.cluster_represent);
}
