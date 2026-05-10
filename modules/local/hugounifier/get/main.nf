process HUGOUNIFIER_GET {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b60674d48ce2e14fec6d6323687b1188e0fe71859f5c82332ea303db69a122d2/data'
        : 'community.wave.seqera.io/library/pip_hugo-unifier:bedd626d591c5003'}"

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
