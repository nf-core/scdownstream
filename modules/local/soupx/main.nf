process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7f7c86baec26b736bdccdc3ee2bcd862b919ced4a82bbe9ef09c646fa6b94117/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-seurat_r-soupx:50aef4aa28333b3c' }"

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
