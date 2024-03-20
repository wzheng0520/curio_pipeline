#!/usr/bin/env python
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pyarrow.compute as pc
import funcy
import sys
import glob
import os
from datetime import datetime

sample_id = sys.argv[1]
r1_databse = sys.argv[2]


def compile_counts(counts, t_counts, tag):
    concatenated = pa.concat_tables([counts, t_counts])
    compiled = concatenated.group_by(tag).aggregate([(tag + "_count", "sum")])
    # compiled = compiled.rename_columns([tag+"_count", tag]) # to be consistent with the count results
    compiled = compiled.rename_columns(
        [tag, tag + "_count"]
    )  # to be consistent with the count results
    # ORDER OF COLUMNS MATTERS
    # BM, BM_count
    return compiled


files = glob.glob(os.path.join(r1_databse, "*.parquet"))
print(len(files))

if len(files) > 2:
    files = files[0:2]
else:
    files = files
print(files)

tag_BB_counts = None
tag_BD_counts = None

for chunk_i, file in enumerate(files):
    print(chunk_i)
    print(file)
    r1_database_i = (
        pq.ParquetDataset(file)
        .read(columns=["r1_proper_structure_matched", "BB", "BD"])
        .combine_chunks()
    )

    # BD proper_structure_matched
    matched_i = r1_database_i.filter((pc.field("r1_proper_structure_matched") == True))
    tag_BB_count_i = matched_i.group_by("BB").aggregate([("BB", "count")])
    tag_BD_count_i = matched_i.group_by("BD").aggregate([("BD", "count")])

    with funcy.print_durations("calculate-BB-counts"):
        if tag_BB_counts is not None:
            tag_BB_counts = compile_counts(tag_BB_counts, tag_BB_count_i, "BB")
        else:
            tag_BB_counts = tag_BB_count_i

    with funcy.print_durations("calculate-BD-counts"):
        if tag_BD_counts is not None:
            tag_BD_counts = compile_counts(tag_BD_counts, tag_BD_count_i, "BD")
        else:
            tag_BD_counts = tag_BD_count_i

tag_BD_counts.to_pandas().to_csv(
    sample_id + ".numProperReads_BD.txt", sep="\t", index=None, header=None
)
pq.write_table(tag_BB_counts, sample_id + ".numProperReads_BB")
