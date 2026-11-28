process DECOUPLER_NETWORK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2821b036063d6024ae937e9f6cbc1b777333ac3bdcca793f68e1fbe868ae976/data'
        : 'community.wave.seqera.io/library/decoupler_activity:374661b94d4383e8'}"

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
