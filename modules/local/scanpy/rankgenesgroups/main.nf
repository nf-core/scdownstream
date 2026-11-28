process SCANPY_RANKGENESGROUPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c75e8f622e621c4de571eb70f474eb4ea1e409723e2cc4488e308d2f32fb36/data'
        : 'community.wave.seqera.io/library/anndata_pandas_pyarrow_python_pruned:c38efc9f3be39e99'}"

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)
    tuple val(filter_col), val(filter_val)
    val(method)
    val(rank_key)
    val(apply_log1p)

    output:
    tuple val(meta), path("*.h5ad")           , emit: h5ad, optional: true
    path "*.pkl"                              , emit: uns, optional: true
    path "*.png"                              , emit: plots, optional: true
    path "*.csv"                              , emit: markers, optional: true
    tuple val(meta), path("*_results.parquet"), emit: results, optional: true
    path "*_mqc.json"                         , emit: multiqc_files, optional: true
    path "versions.yml"                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('rank_genes_groups.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "${prefix}.png"
    touch "${prefix}_dotplot.png"
    touch "${prefix}_mqc.json"
    touch "${prefix}_dotplot_mqc.json"
    touch "${prefix}_results.parquet"
    touch "versions.yml"
    """
}
