#include "bin_map.h"

#include <iostream>
#include <memory>


bin_map::bin_map()
{
    /* init bin_map, calculate bin map start index for each level() */
    _binOffsetsExtended = (binNumType *)calloc(NUM_BIN_LEVELS, sizeof(binNumType));

    // static int binOffsets[] = {4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};
	for (binNumType i=NUM_BIN_LEVELS-2 ; i >= 0; --i) {
		_binOffsetsExtended[i] = _binOffsetsExtended[i+1] + (1 << ((NUM_BIN_LEVELS-2-i) * 3));
	}
}

bin_map::~bin_map() {
	free(_binOffsetsExtended);
}

int bin_map::update_bin_map(const Cluster &query_cluster, const float &max_diff,
    int32_t max_dist, const float &min_overlap)
{
    const std::string &chr = query_cluster.cluster_represent.ref_name1;
    const SVTYPE sv_type = query_cluster.cluster_represent.type;

	main_map_type::iterator mainIter = _main_map.find(chr);
	if (mainIter == _main_map.end()) {
		// given chrom not in map, add query cluster to the bin map
        add_to_bin_map(query_cluster);
		return 1;
	}

    bins_type &type_map = mainIter->second;
    bins_type::iterator typeIter = type_map.find(sv_type);
    if (typeIter == type_map.end()) {
        // given sv type not in type_map
        add_to_bin_map(query_cluster);
        return 1;
    }

    bins &index_map = typeIter->second;

    // expand query cluster 1000bp both side to handle intervals
    // near the bin boundary
    Cluster query_cluster_expand;
    if (query_cluster.cluster_represent.pos1 > 1 << _binFirstShift) {
        query_cluster_expand.cluster_represent.pos1 = query_cluster.cluster_represent.pos1 - 1000;
        query_cluster_expand.cluster_represent.pos2 = query_cluster.cluster_represent.pos2 + 1000;
    // } else if (query_cluster_expand.cluster_represent.pos2 == 0){ // prevent negative end
    //     query_cluster_expand.cluster_represent.pos1 = query_cluster.cluster_represent.pos1;
    //     ++query_cluster_expand.cluster_represent.pos2;
    } else {
        query_cluster_expand.cluster_represent.pos1 = query_cluster.cluster_represent.pos1;
        query_cluster_expand.cluster_represent.pos2 = query_cluster.cluster_represent.pos2;
    } 

    binNumType startPos = (binNumType)(query_cluster_expand.cluster_represent.pos1);
    binNumType endPos = (binNumType)(query_cluster_expand.cluster_represent.pos2);

    binNumType startBin = (startPos >> _binFirstShift);
    binNumType endBin = ((endPos) >> _binFirstShift);


    /* SYNOPSIS:
        1. Loop through each BIN for each SVYTPE of each chrom.
        2. For each BIN, Loop through clusters in the BIN, if the query cluster
           match the cluster in the BIN, update the cluster in the BIN by merge
           it with the query cluster and return 0;
        3. If the query cluster do not match any cluster in all BINs, add it to
           the bin map and return 1;
    */

    for (binNumType i = 0; i < NUM_BIN_LEVELS; ++i) {
        binNumType offset = _binOffsetsExtended[i];
        for (binNumType j = (startBin+offset); j <= (endBin+offset); ++j)  {

        	// move to the next bin if this one is empty
            bins::iterator binsIter = index_map.find(j);
        	if (binsIter == index_map.end()) {
        		continue;
        	}

        	cluster_list &bin = binsIter->second;

        	for (auto iter = bin.begin(); iter != bin.end(); ++iter) {
            	if (query_cluster.intersect(*iter, max_diff, max_dist,
                    min_overlap))
                {
                    // upater cluster in the bin map 
                    // only update the first matched cluster in the bin map
                    iter->merge(query_cluster);
                    // std::cerr << "INFO: merge: " 
                    //     << query_cluster.cluster_represent.id << std::endl;
                    // std::cerr << "INFO: merge: " << iter->ref_name << ";"
                    //     << iter->start << ";"
                    //     << iter->end << ";"
                    //     << iter->type << "\t"
                    //     << query_cluster.ref_name << ";"
                    //     << query_cluster.start << ";"
                    //     << query_cluster.end << ";"
                    //     << query_cluster.type << std::endl;
                    return 0;
            	}
            }
        }
        startBin >>= _binNextShift;
        endBin >>= _binNextShift;
    }
    // new cluster
    add_to_bin_map(query_cluster);
    return 1;
}

void bin_map::get_clusters(std::vector<std::unique_ptr<Cluster>> &clusters) {
    for (auto &map1 : _main_map) {
        for (auto &map2 : map1.second) {
            for (auto &map3 : map2.second) {
                for (auto &_cluster : map3.second) {
                    clusters.push_back(std::unique_ptr<Cluster>(new Cluster(_cluster)));
                }
            }
        }
    }
}

bool bin_map::add_to_bin_map(const Cluster &_cluster)
{
	// Get chr, bin.
	const std::string &chr = _cluster.cluster_represent.ref_name1;
    const SVTYPE &sv_type = _cluster.cluster_represent.type;

    // cast uint32_t to int32_t, make sure that the pos num not exceed 2^31
	binNumType startPos = (binNumType)(_cluster.cluster_represent.pos1);
	binNumType endPos = (binNumType)(_cluster.cluster_represent.pos2);
	binNumType binNum = getBin(startPos, endPos);

	if (binNum < 0 || binNum > NUM_BINS) {
		fprintf(stderr, "ERROR: Received illegal bin number %u from getBin call.\n", binNum);
		return false;
	}

	_main_map[chr][sv_type][binNum].push_back(_cluster);
	return true;
}

bin_map::binNumType bin_map::getBin(binNumType start, binNumType end) const {
    // --end;
    start >>= _binFirstShift;
    end   >>= _binFirstShift;

    for (binNumType i = 0; i < NUM_BIN_LEVELS; ++i) {
        if (start == end) {
        	return _binOffsetsExtended[i] + start;
        }
        start >>= _binNextShift;
        end   >>= _binNextShift;
    }
    //failure
    return -1;
}
