process SCANPY_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scanpy:1.10.1--ea08051addf267ac':
        'community.wave.seqera.io/library/scanpy:1.10.1--0c8c97148fc05558' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    min_genes       = meta.min_genes ?: 1
    min_cells       = meta.min_cells ?: 1
    min_counts_gene = meta.min_counts_gene ?: 1
    min_counts_cell = meta.min_counts_cell ?: 1
    max_mito_fraction   = meta.max_mito_fraction ?: 100
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'filter.py'
}
