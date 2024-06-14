process INTEGRATION_HARMONY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/harmony-pytorch_scanpy:aa1fe21675fae38b':
        'community.wave.seqera.io/library/harmony-pytorch_scanpy:f155fc8e99af0cab' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad") , emit: h5ad
    path "*.pkl"                    , emit: obsm
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'harmony.py'
}
