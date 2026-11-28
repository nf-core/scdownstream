process CELLDEX_FETCHREFERENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f783d952dbf89b36cf16165b825608ee84b0a917cf5250e1515f9bba6d8b665c/data'
        : 'community.wave.seqera.io/library/bioconductor-alabaster.base_bioconductor-alabaster.se_bioconductor-gypsum_bioconductor-hdf5array_pruned:99b54e4c9a06f2d6'}"

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
