process EDGEPYTHON_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbd6f1ae2b21ada0279845e63d73fed5c2c70b874524871d38d8fbad16e815cf/data'
        : 'community.wave.seqera.io/library/edgepython_differential:0be26670466b0c78'}"

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
