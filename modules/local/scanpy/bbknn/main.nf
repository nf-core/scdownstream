process SCANPY_BBKNN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92cfd1aa3ab8886c89cc35fd3d8dd538e92843852f8f4d81f5d406ec4748e9e5/data' :
        'community.wave.seqera.io/library/bbknn_python_pyyaml_scanpy:adc0dbad6ddd2876' }"

    input:
    tuple val(meta), path(h5ad)
    val(batch_col)
    val(input_layer)
    val(log_normalize)

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
