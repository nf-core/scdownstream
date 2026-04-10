process ADATA_READRDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:c4b75a61a89ec006':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:5fae42aabf7a1c5f' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'readrds.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
