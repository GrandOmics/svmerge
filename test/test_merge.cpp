#include "sv.h"
#include "sv_vcf.h"
#include "cluster.h"
#include "bin_map.h"
#include <iostream>
#include <memory>

int main(int argc, char const *argv[])
{
    vcfFile *fp_vcf = vcf_open(argv[1], "r");
    bcf_hdr_t *h_v = bcf_hdr_read(fp_vcf);
    SV _sv;
    bin_map my_bin_map = bin_map();
    int ret;
    int valid = -1;
    // std::cout << _sv.id << std::endl;
    while ((ret = read1_sv_vcf(fp_vcf, h_v, _sv, valid) >=0)) {
        // fprintf(stdout, "%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\n",
        //     _sv.ref_name1.c_str(), _sv.pos1, _sv.ref_name2.c_str(),
        //     _sv.pos2, (int)_sv.type, _sv.length, _sv.sample.c_str(),
        //     _sv.id.c_str());
        if (valid >= 0) {
            Cluster cluster(_sv);
            my_bin_map.update_bin_map(cluster, 0.3, 1000);
        }
    }
    std::vector<Cluster> clusters;
    my_bin_map.get_clusters(clusters);
    for (const auto it: clusters) {
        for (const auto it1: it.SVs) {
            it.cluster_represent.print();
            fprintf(stdout, "\t");
            it1.print();
            fprintf(stdout, "\n");
        }
    }

    bcf_hdr_destroy(h_v);
    vcf_close(fp_vcf);
    return 0;
}