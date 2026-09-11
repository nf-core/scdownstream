process SCANPY_RANKGENESGROUPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a194f12854e4f73b2b4aab27102b3178bfd4dd50fb1332db413205ef3cbfca90/data':
        'community.wave.seqera.io/library/anndata_pandas_pyarrow_python_pruned:f7393f16a921c030' }"

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)
    tuple val(filter_col), val(filter_val)
    val(method)
    val(rank_key)

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
