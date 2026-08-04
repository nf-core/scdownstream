process ADATA_READRDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/28/2863a57333ed9d21fa87f18cb37310d747f777fd4a2b4749298edd4d2d4ed92f/data' :
        'community.wave.seqera.io/library/adata_readrds:771ff9c89e1f71db' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'readrds.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
