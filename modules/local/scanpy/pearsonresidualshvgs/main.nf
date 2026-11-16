process SCANPY_PEARSONRESIDUALS_HVGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0fcad28c1a4f695002836139b6c7c2a65e60794545b7929bc491dbe64b6373e0/data'
        : 'community.wave.seqera.io/library/pyarrow_python_pyyaml_scanpy:bf28df5a8172ceca'}"

    input:
    tuple val(meta), path(h5ad)
    val n_genes
    path excluded_genes
    val batch_key
    val counts_layer

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.parquet", emit: var
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('pearsonresidualshvgs.py')

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
