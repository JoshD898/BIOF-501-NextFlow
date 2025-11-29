#!/usr/bin/env nextflow

params.data_dir = "data"
params.adata_save_file = "adata.h5ad"


workflow {

    download_and_extract_data()
}

process download_and_extract_data {

    output:
    path params.data_dir

    publishDir "${workflow.projectDir}", mode: 'copy'

    script:
    """
    python3 scripts/download_and_extract_data.py --output_dir ${params.data_dir}
    """
}

process create_anndata {

    input:
    path "${params.data_dir}/*"

    output:
    path "${params.adata_save_file}"

    script:
    """
    python3 scripts/create_anndata.py --data_dir ${params.data_dir} --output_file ${params.adata_save_file}
    """
}

