process CELLDEX_FETCHREFERENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/235d42aab096f153f9a10873202d88f2ba8896a3cdb0bac3d8e1d5c2a8fbda6e/data' :
        'community.wave.seqera.io/library/bioconductor-celldex_bioconductor-hdf5array_bioconductor-singlecellexperiment_r-yaml:13bf33457e3e7490' }"

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
