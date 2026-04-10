process CELLDEX_FETCHREFERENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-celldex_bioconductor-hdf5array_bioconductor-singlecellexperiment_r-yaml:c4e76f99d7b45118'
        : 'community.wave.seqera.io/library/bioconductor-celldex_bioconductor-hdf5array_bioconductor-singlecellexperiment_r-yaml:13bf33457e3e7490'}"

    input:
    tuple val(meta), val(ref), val(version)

    output:
    tuple val(meta), path("${prefix}.tar"), emit: tar
    path "versions.yml"                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    template("celldexDownload.R")

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch "${prefix}.tar"
    touch "versions.yml"
    """
}
