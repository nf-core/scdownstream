process HUGOUNIFIER_APPLY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0d8604f4681d37b46850871efd733be8206df7088100a4020f1860a9cea1c1/data'
        : 'community.wave.seqera.io/library/hugo_unifier:2a5a3fafb4d2834f'}"

    input:
    tuple val(meta), path(h5ad, arity: 1), path(changes, arity: 1)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('hugo-unifier'), eval('hugo-unifier --version | grep -oP "(?<=version )[\\d.]+"'), emit: versions_hugo_unifier, topic: versions

    script:
    prefix = task.ext.prefix ?: meta.id
    """
    hugo-unifier apply -i ${h5ad} -c ${changes} -o ${prefix}.h5ad
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.h5ad
    """
}
