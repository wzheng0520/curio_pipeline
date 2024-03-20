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
from Bio import motifs
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
    return compiled


files = glob.glob(os.path.join(r1_databse, "*.parquet"))
print(len(files))

if len(files) > 2:
    files = files[0:2]
else:
    files = files

motifs_proper = None
motifs_improper = None

for chunk_i, file in enumerate(files):
    print(chunk_i)
    print(file)
    r1_database_i = (
        pq.ParquetDataset(file)
        .read(columns=["read1_seq", "proper_structure"])
        .combine_chunks()
    )

    # BD proper_structure_matched
    r1_proper = r1_database_i.filter((pc.field("proper_structure") == True))
    r1_improper = r1_database_i.filter((pc.field("proper_structure") == False))

    t_motifs_proper = pd.DataFrame(
        motifs.create(r1_proper["read1_seq"].to_numpy(), alphabet="ATCGN").counts
    )
    t_motifs_improper = pd.DataFrame(
        motifs.create(r1_improper["read1_seq"].to_numpy(), alphabet="ATCGN").counts
    )

    with funcy.print_durations("calculate-motifs-proper"):
        if motifs_proper is not None:
            motifs_proper = motifs_proper + t_motifs_proper
        else:
            motifs_proper = t_motifs_proper

    with funcy.print_durations("calculate-motifs-improper"):
        if motifs_improper is not None:
            motifs_improper = motifs_improper + t_motifs_improper
        else:
            motifs_improper = t_motifs_improper

positions = list(range(1, motifs_proper.shape[0] + 1))
motifs_proper.insert(loc=0, column="position", value=positions)

positions = list(range(1, motifs_improper.shape[0] + 1))
motifs_improper.insert(loc=0, column="position", value=positions)

motifs_proper.to_csv(sample_id + ".ProperStructure_bybase.txt", index=False, sep="\t")
motifs_improper.to_csv(
    sample_id + ".ImproperStructure_bybase.txt", index=False, sep="\t"
)
