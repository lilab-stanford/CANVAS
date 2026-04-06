#!/usr/bin/env python
# coding: utf-8

import logging
import os
import pandas as pd
import anndata as ad
import scanpy as sc
import scimap as sm
import numpy as np
from concurrent.futures import ProcessPoolExecutor


logging.basicConfig(level=logging.INFO)

logging.info("Importing WS metadata and data")
file_path = '/LUNG.IO_scimap.rawdata.filter.final.csv'
WS_meta = pd.read_csv(
    file_path,
    index_col='patch.id',
    usecols=['patch.id', 'cell.id', 'X_centroid', 'Y_centroid', 'imageid', 'CN_label', 'tu.det_benign'],
    low_memory=False
)

# Import WS data
logging.info("Importing WS data")
WS_count = pd.read_csv(
    '/LUNG.IO_scimap.rawdata.filter.final.csv',
    index_col='patch.id',
    usecols=['patch.id', 'tu.det_benign'],
    low_memory=False
)

adata = ad.AnnData(WS_count)
adata.obs = WS_meta
adata.raw = adata

np.random.seed(42)

adata = sm.tl.spatial_distance(
    adata,
    x_coordinate='X_centroid',
    y_coordinate='Y_centroid',
    phenotype='CN_label',
    subset=None,
    imageid='imageid',
    label='spatial_distance'
)

distance_result = adata.uns['spatial_distance']
distance_output_path = "/dis_results/spatial_distance.csv"
pd.DataFrame(distance_result).to_csv(distance_output_path, index=True)


#
#chmod +x distance_cal.py
#conda activate scimap
#nohup python distance_cal.py > distance_cal.log 2>&1 &

