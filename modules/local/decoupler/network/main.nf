process DECOUPLER_NETWORK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5ca309507e02db5a619c02a9d8bb1d5b969ef6a815cceca6693784ee7fd061ec/data'
        : 'community.wave.seqera.io/library/decoupler_activity:8d8ed42929b0d74b'}"

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
