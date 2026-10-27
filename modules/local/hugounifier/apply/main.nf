process HUGOUNIFIER_APPLY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1fea37c351cbec51dc375e5778be550a5378f842ebdc816c69ff8d7d15c4b8ee/data'
:         'community.wave.seqera.io/library/pip_hugo-unifier:559c1bb3b4aa3398' }"

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
