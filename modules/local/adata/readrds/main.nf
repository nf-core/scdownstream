process ADATA_READRDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/96/96bb395431d12e1db60ae56a0788543a29343cea42657c79d548d015ab69e960/data'
        : 'community.wave.seqera.io/library/adata_readrds:a4ce240ec6f50fd9'}"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('readrds.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
