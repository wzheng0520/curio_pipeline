#!/usr/bin/env python

import dask.dataframe as dd
import pandas as pd
import sys

sample = sys.argv[1]
top_bb_whitelist_parquet = sys.argv[2]

top_bb_whitelist = dd.read_parquet(
    top_bb_whitelist_parquet, columns=["BB", "BM", "BD"]
).compute()

top_bb_whitelist.to_csv(
    sample + ".top_bb_whitelist.txt.gz", index=False, compression="gzip", sep="\t"
)
