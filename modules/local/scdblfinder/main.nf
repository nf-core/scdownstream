process SCDBLFINDER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/709d3cc1733b4f3ddb44a51d7e52dafba6837b495886864bebc53690d95b5f6e/data'
        : 'community.wave.seqera.io/library/scdblfinder:03c14e340018af5c'}"

    input:
    tuple val(meta), path(h5ad), val(dbr), val(batch_col), val(dbr_sd)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.parquet"), emit: predictions
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('scdblfinder.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch versions.yml
    """
}
