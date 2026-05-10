process SEURAT_INTEGRATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b4393c608e642b1232cd7bb84e6c5d7620c4b167462f342a4780307e5e67596b/data':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-rhdf5_r-seurat:71809468c7d8a963' }"

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
