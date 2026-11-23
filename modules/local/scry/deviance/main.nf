process SCRY_DEVIANCE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa26049356fd90f487d0dc0bb8f2491a83e4925696ab6152b3b4c9a87e99c8aa/data'
        : 'community.wave.seqera.io/library/bioconductor-scry_bioconductor-anndatar_bioconductor-singlecellexperiment_bioconductor-rhdf5_r-arrow:2ef29b12964503c9'}"

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
