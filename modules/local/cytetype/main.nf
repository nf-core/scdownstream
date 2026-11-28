process CYTETYPE {
    tag "${meta.id}"
    label 'process_medium'

    // CYTETYPE_API_KEY is attached via conf/modules.config only when params.cytetype_study_context is set.

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/60317d7f1b9e254e4d77eed992c87f833e6ae562ced9b20ad7678805783e6688/data'
        : 'community.wave.seqera.io/library/anndata_pandas_pyarrow_python_pruned:fde96b9090e78ef3'}"

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
    path "${prefix}.parquet", emit: obs
    path "versions.yml", emit: versions, topic: versions

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
