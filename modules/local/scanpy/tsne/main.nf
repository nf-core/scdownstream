process SCANPY_TSNE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04307b777a83d13f004023a046a677adc1c444edcae99a1232cc5b82e88968df/data'
        : 'community.wave.seqera.io/library/leidenalg_pyarrow_python-igraph_python_pruned:085892d5c5a30e57'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val representation

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.parquet", emit: obsm
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_tsne"
    template('tsne.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_tsne"
    """
    touch "${prefix}.h5ad"
    touch "X_${prefix}.parquet"
    touch "versions.yml"
    """
}
