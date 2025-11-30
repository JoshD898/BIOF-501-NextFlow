import scanpy as sc
import scvi
import anndata as ad
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import os

def make_figures(adata: ad.AnnData, model: scvi.model.SCVI, save_dir: str) -> None:
    '''
    Generate and save figures for the given AnnData object using a trained SCVI model.

    Parameters
    ----------
    adata : sc.AnnData
        Input AnnData object containing single-cell data.

    model : scvi.model.SCVI
        Trained SCVI model.

    save_dir : str
        Directory to save the generated figures.
    '''

    adata.obsm["X_scVI"] = model.get_latent_representation()

    # TSNE plot with cell type annotations
    sc.tl.tsne(adata, use_rep="X_scVI")
    sc.settings.figdir = save_dir
    sc.pl.tsne(adata, color="celltypist_labels", title = "Cell Type Annotation", save=f"sne_cell_type.png")

    # Volcano plot for CD14 Monocyte differential expression

    CD14_adata = adata[adata.obs["celltypist_labels"] == "CD14 Monocyte"].copy()

    de = model.differential_expression(
        adata = CD14_adata,
        groupby="disease_status",
        group1="covid_19",
        group2="healthy",
        mode="change",
        delta=0.0,        
    )

    expressed_in_both = de[
    (de["non_zeros_proportion1"] > 0) & 
    (de["non_zeros_proportion2"] > 0)
    ]

    top_genes = expressed_in_both.sort_values("lfc_mean", key=abs, ascending=False).head(10)
    top_gene_names = top_genes.index.tolist()  # SCVI index is gene names

    # Prepare 5x2 figure
    fig, axes = plt.subplots(2, 5, figsize=(20, 12))
    axes = axes.flatten()

    for i, gene in enumerate(top_gene_names):
        ax = axes[i]
        sns.violinplot(
            x="disease_status",
            y=gene,
            data=CD14_adata.to_df()[[gene]].assign(disease_status=CD14_adata.obs["disease_status"]),
            palette={"healthy":"skyblue","covid_19":"salmon"},
            ax=ax
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.set_ylabel("Expression")

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "top10_violin_plots.png"), dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata_file", type=str, help='Path to the AnnData .h5ad file labelled with cell types')
    parser.add_argument("--model_file", type=str, help='Path to the trained scVI model')
    parser.add_argument("--save_dir", type=str, help="Directory to save the generated figures")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.adata_file)
    model = scvi.model.SCVI.load(args.model_file, adata=adata)
    make_figures(adata, model, save_dir=args.save_dir)