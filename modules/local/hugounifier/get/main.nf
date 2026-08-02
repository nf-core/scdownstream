process HUGOUNIFIER_GET {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1fea37c351cbec51dc375e5778be550a5378f842ebdc816c69ff8d7d15c4b8ee/data' :
        'community.wave.seqera.io/library/pip_hugo-unifier:559c1bb3b4aa3398' }"

    input:
    tuple val(meta), val(names), path(h5ads, stageAs: 'input/file_?.h5ad')

    output:
    tuple val(meta), path("${meta.id}/*.csv"), emit: changes
    tuple val("${task.process}"), val('hugo-unifier'), eval('hugo-unifier --version | grep -oP "(?<=version )[\\d.]+"'), emit: versions_hugo_unifier, topic: versions

    script:
    def namedFiles = [[names].flatten(), [h5ads].flatten()].transpose()
    def input = namedFiles.collect { name, h5ad -> "-i ${name}:${h5ad}" }.join(' ')
    """
    hugo-unifier get -o ${meta.id} ${input}
    """

    stub:
    """
    mkdir -p ${meta.id}

    for name in ${names.join(' ')}; do
        touch ${meta.id}/\${name}.csv
    done
    """
}
