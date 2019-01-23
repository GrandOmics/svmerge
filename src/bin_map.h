#ifndef BIN_MAP_H
#define BIN_MAP_H

#include <string>
#include <vector>
#include <map>


#include "sv.h"
#include "cluster.h"


class bin_map {
/*
* using bin index for SV cluster
*/
public:
    bin_map(); // constuctor
    ~bin_map(); // distructor

    /*
    * update bin map by update clusters in it or add new clusters to it
    */
	int update_bin_map(const Cluster &query_cluster,
		const float &max_diff, int32_t max_dist, const float &min_overlap);
    void get_clusters(std::vector<Cluster> &clusters);

private:

    //
    // BIN HANDLING
    //
	typedef int32_t binNumType;
	static const binNumType NUM_BINS = 37449;
	static const binNumType NUM_BIN_LEVELS = 6;

	// bins range in size from 16kb to 512Mb
    // assume that nor chomosome have length exceed 512Mb
	// Bin  0          spans 512Mbp,   # Level 1 # 2^0 bins
	// Bins 1-8        span 64Mbp,     # Level 2 # 2^3 bins
	// Bins 9-72       span 8Mbp,      # Level 3 # 2^6 bins
	// Bins 73-584     span 1Mbp       # Level 4 # 2^9 bins
	// Bins 585-4680   span 128Kbp     # Level 5 # 2^12 bins
	// Bins 4681-37448 span 16Kbp      # Level 6 # 2^15 bins
	binNumType *_binOffsetsExtended;
    // static int binOffsets[] = {4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};
	static const binNumType _binFirstShift = 14;       /* How much to shift to get to finest bin. */
	static const binNumType _binNextShift  = 3;        /* How much to shift to get to next larger bin. */

	typedef std::vector<Cluster> cluster_list; // cluster list
	typedef std::map<binNumType, cluster_list> bins; // for each bin number, have a cluster list
    typedef std::map<SVTYPE, bins> bins_type; // for each type(bit flag) of SV, have a bin map
	typedef std::map<std::string, bins_type> main_map_type; // for each chrom, have bin maps of all SV types.
	main_map_type _main_map; // main data structure to store clusters

	bool add_to_bin_map(const Cluster &_cluster);
	binNumType getBin(binNumType start, binNumType end) const;
	// binNumType getBin(const cluster *_cluster) const;
};

#endif // BIN_MAP_H