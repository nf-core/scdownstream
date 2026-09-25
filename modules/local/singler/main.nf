process CELLTYPES_SINGLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db1479b8f2bde5ec989cb0f7f32f7ae405f16ea3970f4dae91f9ec49de09a7c5/data'
        : 'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-hdf5array_bioconductor-rhdf5_bioconductor-singlecellexperiment_pruned:1a6c1a42a636cbf5'}"

    input:
    tuple val(meta), path(h5ad), val(symbol_col), val(counts_layer)
    tuple val(meta2), val(names), val(labels), path(references)

    output:
    tuple val(meta), path("*_predictions.csv")       , emit: obs
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    tuple val(meta), path("*_distribution.pdf")      , emit: distribution
    tuple val(meta), path("*_heatmap.pdf")           , emit: heatmap
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    counts_layer = counts_layer ?: "X"
    template 'singleR.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_distribution.pdf
    touch ${prefix}_heatmap.pdf
    touch ${prefix}_predictions.csv
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
