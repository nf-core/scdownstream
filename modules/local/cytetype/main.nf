process CYTETYPE {
    tag "${meta.id}"
    label 'process_medium'

    secret secrets.CYTETYPE_API_KEY ? ['CYTETYPE_API_KEY'] : ''
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b1379602a1a1febd0d1715e90506c05aa7945d541538c1f431684be97878dcfa/data'
:         'community.wave.seqera.io/library/python_pyyaml_pydantic_anndata_pruned:87cef5886d0c62d2' }"

    input:
    tuple val(meta), path(h5ad, stageAs: 'input.h5ad')
    val symbol_col
    val study_context
    val group_key
    val rank_key
    val integration
    val resolution

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}.parquet"                   , emit: obs
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
    touch ${prefix}.parquet
    touch versions.yml
    """
}
