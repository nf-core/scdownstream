process ADATA_READRDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:65d6b87e3a6a9270':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:693764013a2a63f5' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path("*.pkl")                  , emit: obsm
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'readrds.py'
}
