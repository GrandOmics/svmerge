#include "sv.h"

SV::SV(std::string _ref_name1, int _pos1, std::string _ref_name2,
    int _pos2, SVTYPE _type, int _length, std::string _sample, std::string _id):
    ref_name1(_ref_name1), pos1(_pos1), ref_name2(_ref_name2),
    pos2(_pos2), type(_type), length(_length), sample(_sample), id(_id)
{
    tra_pos_swap();
}

void SV::tra_pos_swap() {
    if (type == SVTYPE::TRA && ref_name1 > ref_name2) {
        std::swap(ref_name1, ref_name2);
        std::swap(pos1, pos2);
    } 
}

bool SV_intersect(const SV &SV1, const SV &SV2,
    const float &max_diff, int32_t max_dist, const float &min_overlap)
{
    if (SV1.type != SV2.type) {
        return false;
    }

    if (SV1.ref_name1 != SV2.ref_name2 || SV1.ref_name2 != SV2.ref_name2) {
        return false;
    }

    if (SV1.type == SVTYPE::DEL || SV1.type == SVTYPE::DUP ||
        SV1.type == SVTYPE::INV)
    {
        if (SV_overlap_fit(SV1, SV2, min_overlap) &&
            SV_size_fit(SV1, SV2, max_diff) && SV_dist_fit(SV1, SV2, max_dist))
        {
            return true;
        }
    } else if (SV1.type == SVTYPE::TRA) {
        if (SV_dist_fit(SV1, SV2, max_dist)) {
            return true;
        } else {
            fprintf(stderr, "[TRA Merge test] %d.\n", max_dist);
        }
    } else {
        if (SV_size_fit(SV1, SV2, max_diff) && SV_dist_fit(SV1, SV2, max_dist))
        {
            return true;
        }
    }
    return false; 
}

int SV_overlap(const SV &SV1, const SV &SV2)
{
    if (SV1.pos1 > SV2.pos2 || SV1.pos2 < SV2.pos1) {
        return 0;
    } else {
        int overlap_start = SV1.pos1 > SV2.pos1 ? SV1.pos1 : SV2.pos1;
        int overlap_end = SV1.pos2 > SV2.pos2 ? SV2.pos2 : SV1.pos2;
        return overlap_end - overlap_start + 1;
    }
}

bool SV_overlap_fit(const SV &SV1, const SV &SV2, const float &min_overlap)
{
    int overlap = SV_overlap(SV1, SV2);
    if (overlap && (float)overlap/SV1.length >= min_overlap &&
        (float)overlap/SV2.length >= min_overlap)
    {
        return true;    
    } else {
        return false;
    }
}

bool SV_size_fit(const SV &SV1, const SV &SV2, const float max_diff)
{
    if (std::abs(SV1.length - SV2.length)/float(SV1.length) < max_diff &&
        std::abs(SV1.length - SV2.length)/float(SV2.length) < max_diff)
    {
        return true;
    }
    return false;
}

bool SV_dist_fit(const SV &SV1, const SV &SV2, int32_t max_dist)
{
    int32_t dist1 = (int32_t)(SV1.pos1 - SV2.pos1);
    int32_t dist2 = (int32_t)(SV1.pos2 - SV2.pos2);

    if (SV1.length != 0 && 2*SV1.length < max_dist) {
        max_dist = 2*SV1.length;
    }

    if (abs(dist1) < max_dist && abs(dist2) < max_dist) {
        return true;
    }
    return false;
}