process SCANPY_CELLCYCLE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd27aeaf160eaba9a58c029e08f1da74051aa292c2fb043a5dd68fddcde3af93/data'
        : 'community.wave.seqera.io/library/pyyaml_scanpy:3c9e9f631f45553d'}"

    input:
    tuple val(meta), path(h5ad)
    path s_genes
    path g2m_genes
    val symbol_col

    output:
    tuple val(meta), path("${prefix}.pkl"),  emit: obs
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('cellcycle.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pkl
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
