process DECOUPLER_NETWORK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1f600306cb960886d6afe6b911dc74cf521e0521d0cfc16f8a1d4f35711ee2c3/data'
        : 'community.wave.seqera.io/library/decoupler_network:f7895b41fb5156f8'}"

    input:
    tuple val(meta), path(custom_network, stageAs: 'custom/*')
    val species

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: network
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('network.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch versions.yml
    """
}
