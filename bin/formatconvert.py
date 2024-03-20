#!/usr/bin/env python
import pandas as pd
import sys
import os
os.environ[ 'NUMBA_CACHE_DIR' ] = 'tmp'
import scanpy as sc

sample = sys.argv[1]
file_mtx = sys.argv[2]
file_barcodes = sys.argv[3]
file_features = sys.argv[4]
file_matched_whitelist = sys.argv[5]


# rotate tile to match physical orientation
def rotate_whitelist(
    matched_barcode: pd.DataFrame,
) -> pd.DataFrame:
    adjusted_matched_barcode = {
        "SPATIAL_1": -matched_barcode["SPATIAL_2"],
        "SPATIAL_2": matched_barcode["SPATIAL_1"],
    }
    adjusted_matched_barcode = pd.DataFrame(adjusted_matched_barcode)
    return adjusted_matched_barcode


# import expression matrix, barcodes and features
# into an adata object
matched_adata = sc.read_mtx(file_mtx).transpose()
barcodes = pd.read_csv(file_barcodes, header=None)
features = pd.read_csv(file_features, header=None)

matched_adata.var_names = features[0]
matched_adata.obs_names = barcodes[0]

# import matched whitelist
# and rotate to physical tile orientation
matched_bead_whitelist = pd.read_csv(file_matched_whitelist, header=0, index_col=0)
matched_bead_whitelist = rotate_whitelist(matched_bead_whitelist)

# load rotated matched whitelist into the adata object
matched_adata.obsm["spatial"] = matched_bead_whitelist.loc[
    matched_adata.obs_names.tolist(),
].to_numpy()
matched_adata.obsm["X_spatial"] = matched_adata.obsm["spatial"]

# output
matched_adata.write_h5ad(sample + "_anndata.h5ad")
