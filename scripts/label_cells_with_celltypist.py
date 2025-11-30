import celltypist
from celltypist import models
import anndata as ad
import scanpy as sc
import argparse

def annotate_celltypist(adata: ad.AnnData, model_name: str) -> ad.AnnData:
    '''
    Annotate cell types in an AnnData object using CellTypist.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object containing single-cell data.

    model_name : str
        Name of the CellTypist model to use for annotation.
        View all models here: https://www.celltypist.org/models
    '''
    raw_data = adata.X

    sc.pp.normalize_total(adata, target_sum = 10000)
    sc.pp.log1p(adata)

    adata.var_names = adata.var['gene_name'].astype(str).values

    models.download_models(model = model_name)
    predictions = celltypist.annotate(adata, model=model_name, majority_voting=False)

    adata.obs['celltypist_labels'] = predictions.predicted_labels

    adata.X = raw_data

    return adata


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--adata_file", type=str, help="Path to the input AnnData .h5ad file")
    parser.add_argument("--save_path", type=str, help="Path to save the annotated AnnData .h5ad file")
    parser.add_argument("--model_name", type=str, help="CellTypist model name. https://www.celltypist.org/models")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.adata_file)

    adata = annotate_celltypist(adata, model_name=args.model_name)

    adata.write_h5ad(args.save_path)