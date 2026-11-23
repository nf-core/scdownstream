process SCRY_DEVIANCE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45901e03bd44949ad5b11dea802579c7e1d18638086281fcc26202d9db7537f4/data'
        : 'community.wave.seqera.io/library/scry_deviance:6d5e52d3e9cc908f'}"

    input:
    tuple val(meta), path(h5ad)
    val n_genes
    path excluded_genes
    val batch_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*_features.parquet", emit: var
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('deviance.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript -e 'prefix <- gsub(":", "_", "${prefix}", fixed = TRUE); file.create(paste0(prefix, c(".h5ad", "_features.parquet")))'
    touch versions.yml
    """
}
