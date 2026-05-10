process CELLTYPES_SINGLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e1ca6a1732f3fc19cdbcd1ec1872fabb868109dc29fc3e530c1be05d5bc0e5f/data':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-celldex_bioconductor-hdf5array_bioconductor-rhdf5_pruned:74e0a8e51f7ab89c' }"

    input:
    tuple val(meta), path(h5ad), val(symbol_col), val(counts_layer)
    tuple val(meta2), val(names), val(labels), path(references)

    output:
    tuple val(meta), path("*.csv")             , emit: obs
    tuple val(meta), path("*_distribution.pdf"), emit: distribution
    tuple val(meta), path("*_heatmap.pdf")     , emit: heatmap
    path "versions.yml"                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    counts_layer = counts_layer ?: "X"
    template 'singleR.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_distribution.pdf
    touch ${prefix}_heatmap.pdf
    touch ${prefix}.csv
    touch versions.yml
    """
}
