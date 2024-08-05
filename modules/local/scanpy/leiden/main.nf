process SCANPY_LEIDEN {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/leidenalg_python-igraph_scanpy:4c35b8e384d7de2d':
        'community.wave.seqera.io/library/leidenalg_python-igraph_scanpy:ba212a162d290cb6' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.pkl"                   , emit: obs
    path "*.png"                   , emit: plot
    path "*.json"                  , emit: multiqc_files
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    resolution = task.ext.resolution ?: meta.resolution ?: 1.0
    template 'leiden.py'
}
