#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "sv.h"
#include "yutils.h"
#include "htslib/vcf.h"

#include "sv_vcf.h"

BND::BND(const std::string &bnd_str) {
    std::string site_str;
    auto first_right_braket = bnd_str.find('[');
    if (first_right_braket != std::string::npos) {
        if(first_right_braket != 0) {
            type = BND_TYPE::T1;
        } else {
            type = BND_TYPE::T4;
        }
        auto second_right_braket = bnd_str.find('[', first_right_braket+1);
        auto site_str_lenght = second_right_braket - first_right_braket - 1;
        site_str = bnd_str.substr(first_right_braket + 1, site_str_lenght);
    }
    auto first_left_braket = bnd_str.find(']');
    if (first_left_braket != std::string::npos) {
        if (first_left_braket != 0) {
            type = BND_TYPE::T2;
        } else {
            type = BND_TYPE::T3;
        }
        auto second_left_braket = bnd_str.find(']', first_left_braket+1);
        auto site_str_lenght = second_left_braket - first_left_braket - 1;
        site_str = bnd_str.substr(first_left_braket + 1, site_str_lenght);
    }
    std::vector<std::string> site_str_token;
    Tokenize(site_str, site_str_token, ':');
    chrom = site_str_token[0];
    pos = std::stoi(site_str_token[1]);
}

int get_svtype(bcf_hdr_t *h_vcf, bcf1_t *v, std::string &svtype_str) {
    int n1 = 0;
    char *svtype = NULL;
    int ret1 = bcf_get_info_string(h_vcf, v, "SVTYPE", &svtype, &n1);
    if (ret1 < 0) {
        return ret1;
    } 
    svtype_str = svtype;
    free(svtype);
    return 0;
}

int get_sv_end(bcf_hdr_t *h_vcf, bcf1_t *v, int &sv_end) {
    int n1 = 0;
    int32_t *sv_end_ptr = NULL;
    int ret1 = bcf_get_info_int32(h_vcf, v, "END", &sv_end_ptr, &n1);
    if (ret1 < 0) {
        return ret1;
    } 
    sv_end = *sv_end_ptr;
    free(sv_end_ptr);
    return 0;
}

int get_sv_length(bcf_hdr_t *h_vcf, bcf1_t *v, int &sv_length) {
    int n1 = 0;
    int32_t *sv_len_ptr = NULL;
    int ret1 = bcf_get_info_int32(h_vcf, v, "SVLEN", &sv_len_ptr, &n1);
    if (ret1 < 0) {
        return ret1;
    } 
    sv_length = abs(*sv_len_ptr);
    free(sv_len_ptr);
    return 0;
}

// get first sample name in the vcf file
void get_first_vcf_sample(bcf_hdr_t *h_vcf, std::string &sample_name) {
    int nsamples = bcf_hdr_nsamples(h_vcf);
    if (nsamples < 1) {
        fprintf(stderr, "[get_first_vcf_sample] Error! Less than one"
            " sample in vcf header.\n");
        std::exit(1);
    } else if (nsamples > 1) {
        fprintf(stderr, "[get_first_vcf_sample] Warning! More than one"
            " sample in vcf header, only the first one is kept for merge.\n");
    }
    sample_name = h_vcf->samples[0];
}

int read1_sv_vcf(vcfFile *fp_vcf, bcf_hdr_t *h_vcf, SV &_sv, int &valid) {
    bcf1_t *v = bcf_init();

    int ret;
    while((ret = vcf_read1(fp_vcf, h_vcf, v)) >= 0) {
        // sv init
        std::string chrom1;
        int pos1;
        std::string chrom2;

        std::string sample;
        get_first_vcf_sample(h_vcf, sample);

        // chrom and pos
        chrom1 = bcf_hdr_id2name(h_vcf, v->rid);
        pos1 = v->pos + 1; // convert to '1' based coordinate

        // svtype
        int ret1 = 0;
        std::string svtype;
        if ((ret1 = get_svtype(h_vcf, v, svtype)) < 0) {
            fprintf(stderr, "[read_vcf] Error! Can not get SVTYPE from SV: %s,"
                " bcf_get_info_string return %d", v->d.id, ret1);
            ret1 = 0; // reset ret
            std::exit(1);
        }

        int sv_end;
        int sv_length;

        if (svtype == "DEL" || svtype == "DUP" || svtype == "INS"
            || svtype == "INV")
        {
            if ((ret1 = get_sv_end(h_vcf, v, sv_end)) < 0) {
                fprintf(stderr, "[read_vcf] Error! Can not get END from SV: %s,"
                    " bcf_get_info_int32 return %d", v->d.id, ret1);
                ret1 = 0; // reset ret
                std::exit(1);
            }
            if ((ret1 = get_sv_length(h_vcf, v, sv_length)) < 0) {
                fprintf(stderr, "[read_vcf] Error! Can not get SVLEN from SV:"
                " %s, bcf_get_info_int32 return %d", v->d.id, ret1);
                ret1 = 0; // reset ret
                std::exit(1);
            }
            SVTYPE sv_enum_type;
            if (svtype == "DEL") {
                sv_enum_type = SVTYPE::DEL;
            } else if (svtype == "DUP") {
                sv_enum_type = SVTYPE::DUP;
            } else if (svtype == "INS") {
                sv_enum_type = SVTYPE::INS;
            } else {
                sv_enum_type = SVTYPE::INV;
            }
            
            valid = 0;
            _sv = SV(chrom1, pos1, chrom1, sv_end, sv_enum_type,
                sv_length, sample, v->d.id);
        } else if (svtype == "TRA")
        {
            if ((ret1 = get_sv_end(h_vcf, v, sv_end)) < 0) {
                fprintf(stderr, "[read_vcf] Error! Can not get END from SV: %s,"
                    " bcf_get_info_int32 return %d", v->d.id, ret1);
                ret1 = 0; // reset ret
                std::exit(1);
            }
            if ((ret1 = get_sv_length(h_vcf, v, sv_length)) < 0) {
                sv_length = 0; // TRA do not have lenght
            }
            valid = 0;
            _sv = SV(chrom1, pos1, chrom1, sv_end, SVTYPE::TRA,
                sv_length, sample, v->d.id);
        } else if (svtype == "BND")
        {
            // only one allele exist in sv vcf file, 0: ref, 1: allele1
            std::string bnd_str = v->d.allele[1];
            BND bnd(bnd_str);
            sv_end = bnd.pos;
            if (bnd.chrom != chrom1) {
                valid = 0;
                _sv = SV(chrom1, pos1, bnd.chrom, bnd.pos,
                    SVTYPE::TRA, 0, sample, v->d.id);
            } else if (bnd.type == BND_TYPE::T2 || bnd.type == BND_TYPE::T4) {
                if ((ret1 = get_sv_length(h_vcf, v, sv_length)) < 0) {
                    fprintf(stderr, "[read_vcf] Warning! Can not get SVLEN from"
                        " BND(INV) SV: %s, bcf_get_info_int32 return %d",
                        v->d.id, ret1);
                    ret1 = 0; // reset ret
                    sv_length = abs(bnd.pos - pos1);
                }
                valid = 0;
                _sv = SV(chrom1, pos1, bnd.chrom, bnd.pos,
                    SVTYPE::INV, sv_length, sample, v->d.id);
            } else {
                // Drop non-TRA and non-INV BND record
                fprintf(stderr, "[read_vcf] Warning! Drop %s for SV ID:"
                    " %s\n, this BND record is neigther TRA nor INV.",
                    svtype.c_str(), v->d.id);
                valid = -1;
            }
        } else // Drop Complex SV
        {
            fprintf(stderr, "[read_vcf] Warning! Drop complex SV %s for SV ID:"
                " %s\n", svtype.c_str(), v->d.id);
            valid = -1;
        }
        return ret;
    }
    bcf_destroy(v);
    return ret;
}