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
r2_gene_barcode_umi_db = sys.argv[2]


def compile_counts(counts, t_counts, tag):
    try:
        concatenated = pa.concat_tables([counts, t_counts])
    except Exception as e:
        print("Error found")
        print(counts.to_pandas())
        print(t_counts.to_pandas())
        raise e
    compiled = concatenated.group_by(tag).aggregate([(tag + "_count", "sum")])
    compiled = compiled.rename_columns([tag, tag+"_count"])

    #compiled_bm = compiled['BM'].to_numpy()
    #combined_bm_counts = compiled['BM_count_sum'].to_numpy()
    #compiled = pa.Table.from_pydict(
    #    {
    #        'BM': compiled_bm,
    #        'BM_count': combined_bm_counts,
    #    }
    #)
    # compiled = compiled.rename_columns([tag+"_count", tag]) # to be consistent with the count results
    # compiled = compiled.rename_columns(
    #     [f"{tag}_count", tag]
    # )  # to be consistent with the count results

    return compiled


start_time = datetime.now()
print(start_time)

files = glob.glob(os.path.join(r2_gene_barcode_umi_db, "*.parquet"))

total_genic = 0
tag_BM_counts = None

for chunk_i, file in enumerate(files):
    print(chunk_i)
    # umi_db
    r2_gene_barcode_umi_db_i = pq.ParquetDataset(file).read(columns=["BM", "XT"]).combine_chunks()

    #read1_read2_merged_ds = ds.dataset(
    #    file,
    #    format="parquet",
    #    partitioning="hive",
    #)

    #r2_gene_barcode_umi_db_i = (
    #    read1_read2_merged_ds.filter(ds.field("r1_proper_structure_matched") == True)
    #    .filter(ds.field("XT").is_valid())
    #    .filter(ds.field("BM").is_valid())
    #    .to_table(columns=["BM"])
    #)

    total_genic += r2_gene_barcode_umi_db_i.shape[0]

    # BM count for genic reads
    tag_BM_count_i = r2_gene_barcode_umi_db_i.group_by("BM").aggregate([("BM", "count")])

    if tag_BM_counts is not None:
        with funcy.print_durations("calculate-BM-counts"):
            tag_BM_counts = compile_counts(tag_BM_counts, tag_BM_count_i, "BM")
    else:
        tag_BM_counts = tag_BM_count_i

pq.write_table(tag_BM_counts, sample_id + ".numReads_perBM")  ## necessary read2

read2summary = {}
read2summary["Assigned"] = int(total_genic)
read2summary_df = pd.DataFrame({"Values": list(read2summary.keys()), "Counts": list(read2summary.values())})
read2summary_df.to_csv(sample_id + ".read2summary.txt", sep="\t", index=False, header=False)

end_time = datetime.now()
print(end_time)
print(end_time - start_time)
