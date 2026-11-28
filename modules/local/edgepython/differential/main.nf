process EDGEPYTHON_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/af/affdf4c3a8936ecb29ae8472c7ff0e56eef12ccaf37beb820a6ee6f769b1fb07/data'
        : 'community.wave.seqera.io/library/edgepython_sc_differential:83628bf2485ba953'}"

    input:
    tuple val(meta), path(h5ad)
    val reference_condition

    output:
    tuple val(meta), path("${prefix}_*_results.parquet"), emit: results
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stratum_results.parquet
    touch versions.yml
    """
}
