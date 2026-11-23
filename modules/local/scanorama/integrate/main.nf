process SCANORAMA_INTEGRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1f2ae96119e022cecb83d5b3166d783528ce4dd2f999c37b2e76e64daadd6ffa/data'
        : 'community.wave.seqera.io/library/pyarrow_python_pyyaml_scanpy_scanorama:f20e88c532f183ec'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val batch_col
    val input_layer
    val log_normalize

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.parquet", emit: obsm
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('integrate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "X_${prefix}.parquet"
    touch "versions.yml"
    """
}
