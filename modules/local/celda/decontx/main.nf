process CELDA_DECONTX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b65d7d3f33b0423a970eb858e15996b12c901f8277b70ae3feeb7eb3e5e3d868/data':
        'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-celda_bioconductor-rhdf5_bioconductor-singlecellexperiment:ae94f059bb2c3e7f' }"

    input:
    tuple val(meta), path(h5ad), path(raw)
    val(batch_col)
    val(input_layer)
    val(output_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'decontx.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
