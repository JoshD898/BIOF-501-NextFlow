#!/usr/bin/env nextflow

params.data_dir = "data"
params.adata_save_file = "adata.h5ad"
params.scvi_save_file = "scvi_trained_model"
params.max_cells_per_sample = 2000
params.scvi_train_params = "{\"max_epochs\": 30}"
params.celltypist_model = "Adult_COVID19_PBMC.pkl"
params.celltypist_save_file = "celltypist_annotated_adata.h5ad"
params.figures_dir = "figures"
params.output_dir = "output"

workflow {
    data = download_and_extract_data()
    adata = downsample_and_create_anndata(data)
    scvi = train_scvi(adata)
    celltype_adata = label_cells_with_celltypist(adata)
    make_figures(celltype_adata, scvi)
}

process download_and_extract_data {

    output:
    path params.data_dir

    cache 'lenient'

    script:
    """
    python3 ${workflow.projectDir}/scripts/download_and_extract_data.py \
    --output_dir ${params.data_dir}
    """
}

process downsample_and_create_anndata {

    input:
    path data_dir

    cache 'lenient'

    output:
    path params.adata_save_file

    script:
    """
    python3 ${workflow.projectDir}/scripts/downsample_and_create_anndata.py \
        --data_dir ${params.data_dir} \
        --save_path ${params.adata_save_file} \
        --max_cells_per_sample ${params.max_cells_per_sample}
    """
}

process train_scvi {

    input:
    path adata_save_file

    cache 'lenient'

    publishDir "${params.output_dir}", mode: 'copy'

    output:
    path params.scvi_save_file

    script:
    """
    python3 ${workflow.projectDir}/scripts/train_scvi.py \
        --adata_file ${params.adata_save_file} \
        --save_path ${params.scvi_save_file} \
        --train_params '${params.scvi_train_params}'
    """
}

process label_cells_with_celltypist {
    input:
    path adata_save_file

    cache 'lenient'

    publishDir "${params.output_dir}", mode: 'copy'

    output:
    path params.celltypist_save_file

    script:
    """
    python3 ${workflow.projectDir}/scripts/label_cells_with_celltypist.py \
        --adata_file ${params.adata_save_file} \
        --model_name ${params.celltypist_model} \
        --save_path ${params.celltypist_save_file}
    """
}

process make_figures {
    input:
    path celltypist_save_file
    path scvi_save_file

    cache 'lenient'

    output:
    path params.figures_dir

    publishDir "${params.output_dir}", mode: 'copy'

    script:
    """
    python3 ${workflow.projectDir}/scripts/make_figures.py \
        --adata_file ${params.celltypist_save_file} \
        --model_file ${params.scvi_save_file} \
        --save_dir ${params.figures_dir}
    """
}

