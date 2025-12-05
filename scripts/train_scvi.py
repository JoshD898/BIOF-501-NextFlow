import scvi
import anndata as ad
import json
import argparse

def train_scvi(adata: ad.AnnData, **kwargs) -> scvi.model.SCVI:
    '''
    Train an scVI model on the provided AnnData object.

    Parameters
    ----------
    adata : ad.Anndata
        AnnData object containing single-cell RNA-seq data

    save_path : str
        Path to save the trained model

    **kwargs : dict
        Any arguments to pass to model.train(), e.g., max_epochs = 10, early_stopping = True
    '''
    scvi.settings.seed(42)
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    model.train(**kwargs)

    return model

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata_file", type=str, help='Path to the AnnData .h5ad file')
    parser.add_argument("--save_path", type=str, help='Save path for the trained model')
    parser.add_argument("--train_params", type=str, default="{}", help="JSON string of parameters to pass to model.train()")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.adata_file)
    train_params = json.loads(args.train_params)

    model = train_scvi(adata, **train_params)

    model.save(args.save_path)