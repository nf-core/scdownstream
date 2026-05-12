process SCANPY_BBKNN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f065d8b046c4b750e1d302d6d68c181bb779e8bc700ae21f87911565f14178e7/data':
        'community.wave.seqera.io/library/bbknn_pyyaml_scanpy:4cf2984722da607f' }"

    input:
    tuple val(meta), path(h5ad)
    val(batch_col)

    output:
    tuple val(meta), path("*.h5ad") , emit: h5ad
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'bbknn.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
