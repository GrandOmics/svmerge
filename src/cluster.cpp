#include "cluster.h"

Cluster::Cluster(std::shared_ptr<SV> &_sv) {
    cluster_represent = _sv;
    SVs.push_back(_sv);
}

bool Cluster::intersect(const Cluster &other,
    const float max_diff, int32_t max_dist) const {
    if (SV_intersect(*cluster_represent, *other.cluster_represent)) {
        return true;
    }
    return false;
}

void Cluster::merge(const SV &_sv) {
    SVs.push_back(std::make_shared<SV>(_sv));
}

void Cluster::merge(const Cluster &other) {
    SVs.push_back(std::make_shared<SV>(*other.cluster_represent));
}