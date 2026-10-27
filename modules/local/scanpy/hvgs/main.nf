process SCANPY_HVGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf'}"

    input:
    tuple val(meta), path(h5ad)
    val n_hvgs
    path excluded_genes
    val input_layer
    val log_normalize
    val batch_key

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.pkl", emit: var
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
    Path(f"{prefix}.pkl").touch()
    PY
    touch versions.yml
    """
}
