process SCANPY_PCA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4cb85cb4a864e74997f047afb5f16ecd6e106cc8793edfc3e91945501b41bbd/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy_pyarrow:4be15c4d717947be' }"

    input:
    tuple val(meta), path(h5ad)
    val key_added
    val input_layer
    val log_normalize

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.parquet"                 , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_pca"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('pca.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_pca"
    """
    touch ${prefix}.h5ad
    touch X_${prefix}.parquet
    touch versions.yml
    """
}
