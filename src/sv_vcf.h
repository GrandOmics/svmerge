#ifndef SV_VCF_H
#define SV_VCF_H

#include <string>
#include <vector>
#include <memory>

#include "sv.h"
#include "htslib/vcf.h"

enum class BND_TYPE {T1, T2, T3, T4};

class BND {
    public:
    BND() = default;
    BND(const std::string &bnd_str);

    BND_TYPE type;
    std::string chrom;
    int32_t pos;
};

int read1_sv_vcf(vcfFile *fp_vcf, bcf_hdr_t *h_vcf, SV &_sv, int &valid);

#endif // SV_VCF_H