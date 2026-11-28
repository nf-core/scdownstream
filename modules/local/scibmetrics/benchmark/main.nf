process SCIBMETRICS_BENCHMARK {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/94a6bbf6b7f2c1bd7c17a6f6d9c7f195422c6cb160c352d6f08519fc1a4f1f3f/data'
        : 'community.wave.seqera.io/library/python_pyyaml_pip_scib-metrics_umap-learn:a00d5755d869faef'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val max_cells
    val subsample_strategy
    val subsample_seed
    val metric_profile

    output:
    path "${prefix}_metrics.tsv"         , emit: metrics
    path "${prefix}_benchmark_info.json" , emit: benchmark_info
    path "${prefix}_mqc.json"            , emit: multiqc_files
    path "versions.yml"                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('benchmark.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metrics.tsv
    touch ${prefix}_benchmark_info.json
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
