process CELLTYPES_CYTETYPE {
    tag "${meta.id}"
    label 'process_medium'

    secret secrets.CYTETYPE_API_KEY ? ['CYTETYPE_API_KEY'] : ''

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f4ddab9c445762a90a6bc0002f4099b39ae886194a6541cd35929f1c69d26f6/data'
        : 'community.wave.seqera.io/library/leidenalg_python-igraph_python_pyyaml_pruned:dd84aadd1ff0dacd'}"

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
