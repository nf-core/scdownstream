process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07e7961b01e47ed3b939c944e57dc902e635f64b00c32667b1a687d2e996cf30/data':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-seurat_r-soupx:0a27a749423d97ae' }"

    input:
    tuple val(meta), path(h5ad), path(raw)
    val(cluster_resolution)
    val(input_layer)
    val(output_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'soupx.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
