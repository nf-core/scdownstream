process CELLDEX_FETCHREFERENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/342280ffef1e258de52897af412b8123ca62bf9e09198a9dea005258d4ac3da0/data'
:         'community.wave.seqera.io/library/bioconductor-celldex_bioconductor-hdf5array_bioconductor-singlecellexperiment_r-yaml:a79ef050d924368f' }"

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
