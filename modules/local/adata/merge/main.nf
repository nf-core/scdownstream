process ADATA_MERGE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pyyaml_scanpy:158b12038812cf13':
        'community.wave.seqera.io/library/pyyaml_scanpy:61c9ab8e312bbe0a' }"

    input:
    tuple val(meta),  path(h5ads)
    tuple val(meta2), path(base)

    output:
    tuple val(meta), path("*_outer.h5ad")    , emit: outer
    tuple val(meta), path("*_inner.h5ad")    , emit: inner
    tuple val(meta), path("*_integrate.h5ad"), emit: integrate
    path("celltype_predictions/*.csv")       , emit: celltypes
    path "gene_intersection.pkl"             , emit: intersect_genes
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    force_obs_cols = task.ext.force_obs_cols ?: params.force_obs_cols ?: ""
    template 'merge.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_outer.h5ad
    touch ${prefix}_inner.h5ad
    touch ${prefix}_integrate.h5ad
    touch gene_intersection.pkl
    touch versions.yml
    """
}
