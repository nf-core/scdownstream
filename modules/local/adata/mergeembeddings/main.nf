process ADATA_MERGEEMBEDDINGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f8e43bfea40f384db72a635038d6dfa776155850b820ae063a3ab458957b561/data'
        : 'community.wave.seqera.io/library/leidenalg_pyarrow_python-igraph_python_pruned:8631c4d5e6493952'}"

    input:
    tuple val(meta), val(integration_key), path(integrated, stageAs: 'integrated.h5ad'), path(base), path(combined)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.parquet"), emit: obs, optional: true
    path ("X_${prefix}.parquet"), emit: obsm
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('merge_embeddings.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch X_${prefix}.parquet
    touch versions.yml
    """
}
