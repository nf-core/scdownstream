process HUGOUNIFIER_APPLY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f9a59258b22fd6742dce1c44932035843de03273f0466cd4d7a2cec0d9ab107/data' :
        'community.wave.seqera.io/library/pip_hugo-unifier:1545f5b4052d6b2d' }"

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
