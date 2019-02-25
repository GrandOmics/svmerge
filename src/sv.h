#ifndef SV_H
#define SV_H

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <map>
#include <string>

enum class SVTYPE:uint8_t {DEL, INS, DUP, INV, TRA};
const std::map<SVTYPE, std::string>
    svtype_map{{SVTYPE::DEL, "DEL"},
        {SVTYPE::INS, "INS"},
        {SVTYPE::DUP, "DUP"},
        {SVTYPE::INV, "INV"},
        {SVTYPE::TRA, "TRA"}}; 

class SV {
    public:
    SV() = default;
    SV(std::string _ref_name1, int _pos1, std::string _ref_name2,
        int _pos2, SVTYPE _type, int _length, std::string _sample,
        std::string _id, std::string _genotype);

    ~SV(){}

    std::string ref_name1;
    int pos1;
    std::string ref_name2;
    int pos2;
    SVTYPE type;
    int length = 0;

    std::string sample;
    std::string id;

    std::string genotype;

    void tra_pos_swap();

    void print() const{
        fprintf(stderr, "%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s",
            ref_name1.c_str(), pos1, ref_name2.c_str(),
            pos2, svtype_map.find(type)->second.c_str(), length, sample.c_str(),
            id.c_str());
    }
};

bool SV_intersect(const SV &SV1, const SV &SV2, const float &max_diff,
    int32_t max_dist, const float &min_overlap);
int SV_overlap(const SV &SV1, const SV &SV2);
bool SV_overlap_fit(const SV &SV1, const SV &SV2, const float &min_overlap);
bool SV_size_fit(const SV &SV1, const SV &SV2, const float max_diff);
bool SV_dist_fit(const SV &SV1, const SV &SV2, int32_t max_dist);

#endif // SV_H