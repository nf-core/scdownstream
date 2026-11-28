process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/522178976abd5e47ee4e8c22df3b92338d6b2ab3bdc861d2e76959f2819d6dff/data'
        : 'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-seurat_r-soupx:7781aaec9e889d55'}"

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
