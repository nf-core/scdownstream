process CELLTYPES_CYTETYPE {
    tag "${meta.id}"
    label 'process_medium'

    secret 'CYTETYPE_API_KEY'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1da96b02584d1ff84d0f4a00be1bece53aa52baacad2f6dbca0bab5584203aa2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy_pip_cytetype:32269855682fb82e'}"

    input:
    tuple val(meta), path(h5ad), val(symbol_col)
    val study_context
    val leiden_resolution
    val n_top_genes

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.pkl") , emit: obs
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('cytetype.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    touch versions.yml
    """
}
