#!/usr/bin/env python

import dask.dataframe as dd
import pandas as pd
import sys

sample = sys.argv[1]
read1_db = sys.argv[2]
top_bb_whitelist_merged_parquet = sys.argv[3]

R1_structure_tags = dd.read_parquet(
    read1_db,
    columns=["read_id", "proper_structure", "cellular_barcode", "keep_illumina"],
)
BB_vs_BM = dd.read_parquet(
    top_bb_whitelist_merged_parquet,
    columns=["cellular_barcode", "molecular_barcode", "BB", "BM", "BD"],
)

R1_structure_tags_merged = dd.merge(
    R1_structure_tags, BB_vs_BM, on=["cellular_barcode"], how="left"
)
R1_structure_tags_merged_persisted = R1_structure_tags_merged.persist()
R1_properstructure_matched = R1_structure_tags_merged_persisted[
    (R1_structure_tags_merged_persisted.BD >= 0)
    & (R1_structure_tags_merged_persisted.proper_structure == True)
]

reads_summary = {}
reads_summary["total"] = R1_structure_tags_merged_persisted.shape[0].compute()
reads_summary["proper_structure"] = (
    R1_structure_tags_merged_persisted["proper_structure"].sum().compute()
)
reads_summary["keep_illumina"] = (
    R1_structure_tags_merged_persisted["keep_illumina"].sum().compute()
)
reads_summary["properstructure_matched"] = R1_properstructure_matched.shape[0].compute()

reads_summary_df = pd.DataFrame(
    {"Values": list(reads_summary.keys()), "Counts": list(reads_summary.values())}
)

numProperReads_BD = R1_properstructure_matched.BD.value_counts().compute()
numProperReads_BM = R1_properstructure_matched.BM.value_counts().compute()
numProperReads_BB = R1_properstructure_matched.BB.value_counts().compute()

reads_summary_df.to_csv(
    sample + ".read_summary.txt", sep="\t", index=False, header=False
)
numProperReads_BD.to_csv(
    sample + ".numProperReads_BD.txt", sep="\t", index=True, header=False
)
numProperReads_BM.to_csv(
    sample + ".numProperReads_BM.txt", sep="\t", index=True, header=False
)
numProperReads_BB.to_csv(
    sample + ".numProperReads_BB.txt", sep="\t", index=True, header=False
)
