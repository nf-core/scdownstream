process SCANPY_LEIDEN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/igraph_scanpy:b5cb16ff6177c14e':
        'community.wave.seqera.io/library/igraph_scanpy:df2be824a2c1834c' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.pkl"                   , emit: obs
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    resolution = task.ext.resolution ?: meta.resolution ?: 1.0
    template 'leiden.py'
}
