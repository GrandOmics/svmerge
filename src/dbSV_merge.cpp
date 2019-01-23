#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>

#include "sv.h"
#include "sv_vcf.h"
#include "cluster.h"
#include "bin_map.h"

void usage() {
    std::cerr << "dbSV merge tool. This tool is used for merge SVs from different sample(callset).\n" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "-h, --help                    print this message and exit\n"
        << "-f, --fofn, FILE              input file, one vcf file path per line [default: None]\n"
        << "-o, --output, FILE            output file [default: None]\n"
        << "-d, --max_distance, INT       maximum distance on both start and end posion of the SV [default: 1000]\n"
        << "-l, --max_length_diff, FLOAT  maximum SV length difference. [default 0.5]\n"
        << "-r, --min_overlap, FLOAT      minimum SV overlap. [default 0.5]\n"
        << "-V, --version                 print version"
        << std::endl;
}

void sv_merge_cluster1(const char *vcf_fn, const float &max_diff,
    int max_dist, const float &min_overlap, bin_map &_bin_map)
{
    vcfFile *fp_vcf = vcf_open(vcf_fn, "r");
    bcf_hdr_t *h_v = bcf_hdr_read(fp_vcf);
    SV _sv;
    int ret;
    int valid = -1;
    while ((ret = read1_sv_vcf(fp_vcf, h_v, _sv, valid) >=0)) {
        if (valid >= 0) {
            Cluster cluster(_sv);
            _bin_map.update_bin_map(cluster, max_diff, max_dist, min_overlap);
        }
    }
    bcf_hdr_destroy(h_v);
    vcf_close(fp_vcf);
}

void sv_merge(const char *vcf_fofn, char *out_fn, const float &max_diff,
    int max_dist, const float &min_overlap)
{
    std::vector<std::string> vcfs;
    std::ifstream vcf_fofn_fp;
    vcf_fofn_fp.open(vcf_fofn);
    std::string vcf_fn;
    while (getline(vcf_fofn_fp, vcf_fn)) {
        vcfs.push_back(vcf_fn);
    }

    bin_map _bin_map = bin_map();
    
    for (const auto it: vcfs) {
        sv_merge_cluster1(it.c_str(), max_diff, max_dist, min_overlap,
            _bin_map);
    }

    FILE *out_fp = fopen(out_fn, "w");
    std::vector<std::unique_ptr<Cluster>> clusters;
    _bin_map.get_clusters(clusters);
    int n = 0;
    for (std::vector<std::unique_ptr<Cluster>>::size_type i = 0; i < clusters.size(); ++i) {
        ++n;
        for (const auto it: clusters[i]->SVs) {
            fprintf(out_fp, "dbsv%d\t", n);
            clusters[i]->cluster_represent.print(out_fp);
            fprintf(out_fp, "\t");
            it.print(out_fp);
            fprintf(out_fp, "\n");
        }
    }
    fclose(out_fp);
}

int main(int argc, char *argv[])
{
    if (argc <= 1) {
        usage();
        return 0;
    }

    std::string __version__ = "v0.0.1";

    // long_option array
    static const struct option long_options[] = {
        {"fofn", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
        {"max_distance", required_argument, 0, 'd'},
        {"max_length_diff", required_argument, 0, 'l'},
        {"min_overlap", required_argument, 0, 'r'},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'V'}
    };

    int c, long_idx;
    const char *opt_str = "f:o:d:l:hV";

    char *vcf_fofn;
    char *output_fn;

    // defaults
    int max_distance = 1000;
    float max_length_diff = 0.5;
    float max_overlap = 0.5;

    // int getopt_long (int argc, char *const *argv, const char *shortopts,
    // const struct option *longopts, int *indexptr)
    while ((c = getopt_long(argc, argv, opt_str, long_options,
        &long_idx)) != -1)
    {
        if (c == 'f') {
            vcf_fofn = optarg;
        } else if (c == 'o') {
            output_fn = optarg;
        } else if (c == 'd') {
            max_distance = atoi(optarg);
        } else if (c == 'l') {
            max_length_diff = atof(optarg);
        } else if (c == 'r') {
            max_overlap = atof(optarg);
        } else if (c == 'h') {
            usage();
            return 0;
        }
        else if (c == 'V') {
            std::cout << "dbSV merge version: " << __version__ << std::endl;
            return 0;
        }
        else {
            std::cerr << "Invalid Arguments." << std::endl;
            std::exit(1);
        }
    }
    sv_merge(vcf_fofn, output_fn, max_length_diff, max_distance, max_overlap);
}