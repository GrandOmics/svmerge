# GrandOmics SV merge tool.

SV merge tool for build dbSV.

## Dependency

**Linux**

cmake3.2+

gcc4.8+

htslib

## Install

```
git clone --recursive git@gitlab.com:archieyoung/svmerge.git
cd svmerge
mkdir build
cd build
cmake ..
make
```

A excutable file `dbSV_merge` will be generated in build/bin directory. You may copy `dbSV_merge` to your `PATH`.

## Usage

```
dbSV_merge

dbSV merge tool. This tool is used for merge SVs from different sample(callset).

Options:
-h, --help                    print this message and exit
-f, --fofn, FILE              input file, one vcf file path per line [default: None]
-o, --output, FILE            output file [default: None]
-d, --max_distance, INT       maximum distance on both start and end posion of the SV [default: 1000]
-l, --max_length_diff, FLOAT  maximum SV length difference. [default 0.5]
-r, --min_overlap, FLOAT      minimum SV overlap. [default 0.5]
-V, --version                 print version
```

## Output

Output file is a simple table seperated by *tab*. Fields are:

ID: merged sv id. start with *dbsv*, end with the rank of the merged sv in the table.

sv fields: CHROM1, POS1, CHROM2, POS2, SVTYPE, SVLEN, SAMPLE_ID, SV_ID

For example:

```
dbsv1   1       725164  1       224200278       DEL     223475114       SAMPLE1       130
dbsv1   1       725254  1       224200328       DEL     223475074       SAMPLE2       63
dbsv1   1       725350  1       224200387       DEL     223475037       SAMPLE3       109
dbsv1   1       725427  1       224200279       DEL     223474852       SAMPLE4       62
dbsv1   1       725268  1       224200315       DEL     223475047       SAMPLE5       107
dbsv1   1       725261  1       224200342       DEL     223475081       SAMPLE6       113
dbsv1   1       725164  1       224200150       DEL     223474986       SAMPLE7       87
dbsv1   1       725263  1       224200280       DEL     223475017       SAMPLE8       86
dbsv1   1       724882  1       224200339       DEL     223475457       SAMPLE9       123
dbsv1   1       725269  1       224200333       DEL     223475064       SAMPLE10      84
```

The first line of one merged sv is the `represent` SV for the cluster.

## Future Plan

Update merged SV table progressively.