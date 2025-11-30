import os
import re
import scanpy as sc
from scipy.io import mmread
import anndata as ad
import pandas as pd
import numpy as np
import argparse
from typing import Optional

def load_anndata(data_dir: str, max_cells_per_sample: Optional[int] = None) -> ad.AnnData:
    '''
    Load the GSE166992 dataset from the data dir, perform downsampling, and create an AnnData object.

    The max_cells_per_sample parameter can be used to limit the number of cells loaded per sample for testing purposes.

    Cells are labelled as 'healthy' or 'covid_19' based on hard-coded sample IDs.

    Healthy sample IDs were taken from GEO: 
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166992

    Parameters
    ----------
    data_dir : str
        Directory containing extracted data files

    max_cells_per_sample : Optional[int]
        Optional parameter to limit the number of cells to load per sample
    
    '''
    HEALTHY_SAMPLES = ['GSM5090454', 'GSM5090448', 'GSM5090446']

    samples = {}

    # fancy regex ensures this works even if files are in random order
    for file in os.listdir(data_dir):
        match = re.match(r"(GSM\d+_[^_]+_5DGE)_(barcodes|features|matrix)\.(tsv|mtx)", file)
        sample_id = match.group(1)
        file_type = match.group(2)
        samples.setdefault(sample_id, {})[file_type] = os.path.join(data_dir, file)

    adatas = []

    for sample, data_paths in samples.items():

        matrix = mmread(data_paths["matrix"]).tocsr().T

        features = pd.read_csv(data_paths["features"], sep="\t", header=None).astype(str)
        features.columns = ['gene_id', 'gene_name', 'feature_type']
        features.index = features["gene_id"].astype(str)

        adata = ad.AnnData(
            X=matrix,
            var = features
        )

        sample_id = sample[:10]
        adata.obs["sample_id"] = sample_id
        
        if sample_id in HEALTHY_SAMPLES:
            adata.obs["disease_status"] = "healthy"
        else:
            adata.obs["disease_status"] = "covid_19"


        ## DOWNSAMPLING
        if max_cells_per_sample is not None:
            np.random.seed(42)
            n = min(max_cells_per_sample, adata.n_obs)
            sampled_indices = np.random.choice(adata.n_obs, size=n, replace=False)
            adata = adata[sampled_indices].copy()

        adatas.append(adata)

    adata = ad.concat(adatas, join='outer', axis=0)
    adata.var = adatas[0].var

    return adata

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, help='Directory containing extracted data files')
    parser.add_argument('--save_path', type=str, help='Save path for AnnData object')
    parser.add_argument('--max_cells_per_sample', type=int, default=None, help='Maximum number of cells to load per sample')
    args = parser.parse_args()

    adata = load_anndata(args.data_dir, args.max_cells_per_sample)
    adata.write(args.save_path)