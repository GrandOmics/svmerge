#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <stdio.h>

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

int get_chr2(bcf_hdr_t *h_vcf, bcf1_t *v, std::string &chr2_str) {
    int n1 = 0;
    char *chr2 = NULL;
    int ret1 = bcf_get_info_string(h_vcf, v, "CHR2", &chr2, &n1);
    if (ret1 < 0) {
        return ret1;
    } 
    chr2_str = chr2;
    free(chr2);
    return 0;
}

int get_sv_end(bcf_hdr_t *h_vcf, bcf1_t *v, int &sv_end) {
    int n1 = 0;
    int32_t *sv_end_prt = NULL;
    int ret1 = bcf_get_info_int32(h_vcf, v, "END", &sv_end_prt, &n1);
    if (ret1 < 0) {
        return ret1;
    } 
    sv_end = *sv_end_prt;
    free(sv_end_prt);
    return 0;
}

int get_sv_length(bcf_hdr_t *h_vcf, bcf1_t *v, int &sv_length) {
    int n1 = 0;
    int32_t *sv_len_prt = NULL;
    int ret1 = bcf_get_info_int32(h_vcf, v, "SVLEN", &sv_len_prt, &n1);
    if (ret1 < 0) {
        return ret1;
    } 
    sv_length = abs(*sv_len_prt);
    free(sv_len_prt);
    return 0;
}

int get_sv_genotype(bcf_hdr_t *h_vcf, bcf1_t *v, std::string &sv_gt) {
    
    int n1 = 0;
    int32_t *gt_arr = NULL;
    
    int ret1 = bcf_get_genotypes(h_vcf, v, &gt_arr, &n1);
    
    if (ret1 < 0) {
        return ret1;
    }
    int nsmpl = bcf_hdr_nsamples(h_vcf);

    int max_ploidy = ret1/nsmpl;
    // gt of the first sample
    int32_t *ptr = gt_arr;
    char aChar[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    int is_phased = 0;
    std::vector<char> gt_tmp;
    for (int i = 0; i < max_ploidy; ++i) {
        if (ptr[i]==bcf_int32_vector_end) {
            break;
        }
        if (bcf_gt_is_missing(ptr[i])) {
            gt_tmp.push_back('.');
            continue;
        }
        int allele_index = bcf_gt_allele(ptr[i]);

        is_phased = bcf_gt_is_phased(ptr[i]);

        gt_tmp.push_back(aChar[allele_index]);
    }

    char sep = '/';
    if (is_phased) {
        sep = '|';
    }

    std::string gt_tmp_str;
    for (auto begin = gt_tmp.begin(); begin != gt_tmp.end() - 1; ++begin) {
        gt_tmp_str += *begin;
        gt_tmp_str += sep;
    }
    gt_tmp_str += *(gt_tmp.end()-1);

    sv_gt = gt_tmp_str;

    free(gt_arr);
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
            " sample in vcf header, only the first one will be kept for merge.\n");
    }
    sample_name = h_vcf->samples[0];
}

int read1_sv_vcf(vcfFile *fp_vcf, bcf_hdr_t *h_vcf, SV &_sv, int &valid) {

    valid = -1; // default is not valid

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

        int ret1 = 0;
        // svtype
        std::string svtype;
        if ((ret1 = get_svtype(h_vcf, v, svtype)) < 0) {
            fprintf(stderr, "[get_svtype] Error! Can not get SVTYPE from SV: "
                "%s, bcf_get_info_string return %d.\n", v->d.id, ret1);
            ret1 = 0; // reset ret
            bcf_destroy(v);
            std::exit(1);
        }
        
        int sv_end;
        int sv_length;

        std::string _genotype;
        if ((ret1 = get_sv_genotype(h_vcf, v, _genotype) < 0)) {
            fprintf(stderr, "[read1_sv_vcf] Warning! Can not get GT from "
                "SV: %s, bcf_get_genotypes return %d.\n", v->d.id, ret1);
            ret1 = 0; // reset ret
            _genotype = "./."; // set gt as "./."
        }

        if (svtype == "DEL" || svtype == "DUP" || svtype == "INS"
            || svtype == "INV")
        {
            if ((ret1 = get_sv_end(h_vcf, v, sv_end)) < 0) {
                fprintf(stderr, "[read1_sv_vcf] Error! Can not get END from SV:"
                    " %s, bcf_get_info_int32 return %d.\n", v->d.id, ret1);
                ret1 = 0; // reset ret
                bcf_destroy(v);
                std::exit(1);
            }
            if ((ret1 = get_sv_length(h_vcf, v, sv_length)) < 0) {
                fprintf(stderr, "[read1_sv_vcf] Error! Can not get SVLEN from "
                    "SV: %s, bcf_get_info_int32 return %d.\n", v->d.id, ret1);
                ret1 = 0; // reset ret
                bcf_destroy(v);
                std::exit(1);
            }
            SVTYPE sv_enum_type;
            if (svtype == "DEL") {
                sv_enum_type = SVTYPE::DEL;
            } else if (svtype == "DUP") {
                sv_enum_type = SVTYPE::DUP;
            } else if (svtype == "INS") {
                sv_enum_type = SVTYPE::INS;
                sv_end = pos1; // force INS END = INS START to fix sniffles INS END bug
            } else {
                sv_enum_type = SVTYPE::INV;
            }
            
            _sv = SV(chrom1, pos1, chrom1, sv_end, sv_enum_type,
                sv_length, sample, v->d.id, _genotype);
            valid = 0;
        } else if (svtype == "TRA")
        {
            if ((ret1 = get_sv_end(h_vcf, v, sv_end)) < 0) {
                fprintf(stderr, "[read1_sv_vcf] Error! Can not get END from SV:"
                    " %s, bcf_get_info_int32 return %d.\n", v->d.id, ret1);
                ret1 = 0; // reset ret
                bcf_destroy(v);
                std::exit(1);
            }
            if ((ret1 = get_sv_length(h_vcf, v, sv_length)) < 0) {
                sv_length = -1; // TRA do not have lenght
                ret1 = 0; // reset ret
            }

            if ((ret1 = get_chr2(h_vcf, v, chrom2)) < 0) {
                fprintf(stderr, "[read1_sv_vcf] Error! Can not get CHR2 from SV:"
                    " %s, bcf_get_info_string return %d.\n", v->d.id, ret1);
                ret1 = 0; // reset ret
                bcf_destroy(v);
                std::exit(1);
            }

            _sv = SV(chrom1, pos1, chrom2, sv_end, SVTYPE::TRA,
                sv_length, sample, v->d.id, _genotype);
            valid = 0;
        } else if (svtype == "BND")
        {
            // only one allele exist in sv vcf file, 0: ref, 1: allele1
            std::string bnd_str = v->d.allele[1];
            BND bnd(bnd_str);
            sv_end = bnd.pos;
            if (bnd.chrom != chrom1) {
                _sv = SV(chrom1, pos1, bnd.chrom, bnd.pos,
                    SVTYPE::TRA, 0, sample, v->d.id, _genotype);
                valid = 0;
            } else if (bnd.type == BND_TYPE::T2 || bnd.type == BND_TYPE::T4) {
                if ((ret1 = get_sv_length(h_vcf, v, sv_length)) < 0) {
                    // fprintf(stderr, "[read1_sv_vcf] Warning! Can not get SVLEN"
                    //     " from BND(INV) SV: %s, bcf_get_info_int32 return %d.\n",
                    //     v->d.id, ret1);
                    ret1 = 0; // reset ret
                    sv_length = abs(bnd.pos - pos1);
                }
                _sv = SV(chrom1, pos1, bnd.chrom, bnd.pos,
                    SVTYPE::INV, sv_length, sample, v->d.id, _genotype);
                valid = 0;
            } else {
                // Drop non-TRA and non-INV BND record
                fprintf(stderr, "[read1_sv_vcf] Warning! Drop %s for SV ID:"
                    " %s\n, this BND record is neigther TRA nor INV.\n",
                    svtype.c_str(), v->d.id);
            }
        }
        bcf_destroy(v);
        return ret;
    }
    bcf_destroy(v);
    return ret;
}