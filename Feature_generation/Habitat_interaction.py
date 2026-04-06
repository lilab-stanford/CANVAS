#!/usr/bin/env python
# coding: utf-8

import logging
import os
import pandas as pd
import anndata as ad
import scimap as sm
import numpy as np
from concurrent.futures import ProcessPoolExecutor

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')

def load_data():
    file_path = '/LUNG.IO_scimap.rawdata.filter.final.csv'
    WS_meta = pd.read_csv(
        file_path,
        index_col='patch.id',
        usecols=['patch.id', 'cell.id', 'X_centroid', 'Y_centroid', 'imageid', 'CN_label', 'tu.det_benign'],
        low_memory=False
    )
    WS_count = pd.read_csv(
        file_path,
        index_col='patch.id',
        usecols=['patch.id', 'tu.det_benign'],
        low_memory=False
    )
    adata = ad.AnnData(WS_count)
    adata.obs = WS_meta
    adata.raw = adata
    return adata

def calculate_spatial_interaction(radius):
    try:
        logging.info(f"Processing radius: {radius}")
        adata = load_data()
        adata_temp = sm.tl.spatial_interaction(
            adata,
            method='radius',
            radius=radius,
            pval_method='zscore',
            permutation=1000,
            phenotype='CN_label',
            imageid='imageid',
            x_coordinate='X_centroid',
            y_coordinate='Y_centroid',
            label='spatial_interaction_radius'
        )
        result_path = f"/results/spatial_interaction_radius_{radius}.csv"
        pd.DataFrame(adata_temp.uns['spatial_interaction_radius']).to_csv(result_path)
        logging.info(f"Completed radius: {radius}")
    except Exception as e:
        logging.error(f"Failed to process radius {radius} with error {e}")

np.random.seed(42)

radii = range(100, 1001, 100)
with ProcessPoolExecutor(max_workers=30) as executor:
    executor.map(calculate_spatial_interaction, radii)


#chmod +x CN_interaction.py
#conda activate scimap
#nohup python CN_interaction.py > CN_interaction.log 2>&1 &

