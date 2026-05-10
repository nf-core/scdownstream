process HUGOUNIFIER_APPLY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b60674d48ce2e14fec6d6323687b1188e0fe71859f5c82332ea303db69a122d2/data'
        : 'community.wave.seqera.io/library/pip_hugo-unifier:bedd626d591c5003'}"

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
