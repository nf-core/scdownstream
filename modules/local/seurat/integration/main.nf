process SEURAT_INTEGRATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-singlecellexperiment_r-seurat:4e7698b0ad6ef366':
        'community.wave.seqera.io/library/bioconductor-singlecellexperiment_r-seurat:fa62c9e0fa2cbc0c' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'integration.R'
}
