process ADATA_TORDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c002326219227ab9a82df4eb816e1883417716ec2c51193c77ecf4ecb8e967f2/data'
        : 'community.wave.seqera.io/library/adata_tords:221549603ba837e4'}"

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
