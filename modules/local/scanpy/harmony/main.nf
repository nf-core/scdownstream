process SCANPY_HARMONY {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ ext.use_gpu ? 'docker.io/nicotru/rapids-singlecell:0.0.1' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/harmonypy_scanpy:25dd7b6ef4c76098':
        'community.wave.seqera.io/library/harmonypy_scanpy:e9c9a621297da9ea' }"

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
