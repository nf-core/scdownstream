process SCANPY_HVGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04307b777a83d13f004023a046a677adc1c444edcae99a1232cc5b82e88968df/data'
        : 'community.wave.seqera.io/library/leidenalg_pyarrow_python-igraph_python_pruned:085892d5c5a30e57'}"

    input:
    tuple val(meta), path(h5ad)
    val n_hvgs
    path excluded_genes
    val input_layer
    val log_normalize
    val batch_key

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.parquet", emit: var
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('hvgs.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<'PY'
    from pathlib import Path
    prefix = "${prefix}".replace(":", "_")
    Path(f"{prefix}.h5ad").touch()
    Path(f"{prefix}.parquet").touch()
    PY
    touch versions.yml
    """
}
