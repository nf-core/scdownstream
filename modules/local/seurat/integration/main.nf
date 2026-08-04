process SEURAT_INTEGRATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba95c9e79901d212e616821dccbabd8998d5620e7bd9cf2ab804df0f96bf037b/data' :
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-glmgampoi_bioconductor-rhdf5_r-seurat:219408537f886c16' }"

    input:
    tuple val(meta), path(h5ad)
    val(batch_col)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('integration.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
