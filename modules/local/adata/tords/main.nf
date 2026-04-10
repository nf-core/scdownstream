process ADATA_TORDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e2d0625efb46bae98303a0e685efa0ee259b939a6490eef13218c16c601bc9b/data':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_bioconductor-singlecellexperiment:b7b9571d025f377e' }"

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
