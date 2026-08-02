process ADATA_TORDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63ce18f987c08fe48ed9a5ec7fee64561a63297771bda7bcb6d9f591162be4c3/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_bioconductor-singlecellexperiment:7f4fd121bcbeb705' }"

    input:
    tuple val(meta), path(h5ad)
    val counts_layer

    output:
    tuple val(meta), path("*.rds"), emit: rds
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    counts_layer = counts_layer ?: 'X'
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'tords.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds
    touch versions.yml
    """
}
