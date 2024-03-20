#!/usr/bin/env python

import pysam
import pandas as pd
import sys


def count_tag(
    read: pysam.libcalignedsegment.AlignedSegment,
    tag: str,
    tag_dict: dict,
):
    tag_value = read.get_tag(tag)
    if tag_value in tag_dict:
        tag_dict[tag_value] += 1
    else:
        tag_dict[tag_value] = 1


def convert_dict_to_df(
    tag_dict: dict,
):
    tag_df = pd.DataFrame(
        {"Values": list(tag_dict.keys()), "Counts": list(tag_dict.values())}
    )
    tag_df = tag_df.sort_values("Counts", ascending=False)
    return tag_df


sample = sys.argv[1]
proper_structure_matched_sorted_bam = sys.argv[2]

bam_proper_structure_matched = pysam.AlignmentFile(
    proper_structure_matched_sorted_bam, "rb"
)
bam_proper_structure_matched_parsed = bam_proper_structure_matched.fetch

BM_dict = {}
BB_dict = {}
BD_dict = {}
XS_dict = {}

for read in bam_proper_structure_matched.fetch():
    if read.mapping_quality > 10:
        count_tag(read, "BM", BM_dict)
        count_tag(read, "BB", BB_dict)
        count_tag(read, "BD", BD_dict)
        count_tag(read, "XS", XS_dict)

BM_df = convert_dict_to_df(BM_dict)
BB_df = convert_dict_to_df(BB_dict)
BD_df = convert_dict_to_df(BD_dict)
XS_df = convert_dict_to_df(XS_dict)

BM_df.to_csv(
    sample + ".numReads_BM_mq_10.matched.txt", sep="\t", index=False, header=False
)
BB_df.to_csv(
    sample + ".numReads_BB_mq_10.matched.txt", sep="\t", index=False, header=False
)
BD_df.to_csv(
    sample + ".numReads_BD_mq_10.matched.txt", sep="\t", index=False, header=False
)
XS_df.to_csv(
    sample + ".numReads_XS_mq_10.matched.txt", sep="\t", index=False, header=False
)
