process CUSTOM_AGGREGATE_PER_CELL_ANNOTATIONS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4cb85cb4a864e74997f047afb5f16ecd6e106cc8793edfc3e91945501b41bbd/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy_pyarrow:4be15c4d717947be' }"

    input:
    tuple val(meta), path(h5ad)
    val cluster_col
    val annotation_col
    val integration
    val resolution

    output:
    tuple val(meta), path("${prefix}.h5ad")            , emit: h5ad
    tuple val(meta), path("${prefix}.parquet")             , emit: obs
    path "${prefix}_aggregate_per_cell_annotations.csv", emit: aggregate_per_cell_annotations, optional: true
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('aggregate_per_cell_annotations.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch versions.yml
    """
}
