process SCDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-scds:1.18.0--b0910f04d88fb193':
        'community.wave.seqera.io/library/bioconductor-scds:1.18.0--aaf652129cf65197' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.csv"), emit: predictions
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'scds.R'
}
