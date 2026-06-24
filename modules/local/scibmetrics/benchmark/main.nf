process SCIBMETRICS_BENCHMARK {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6bd40433b5f1000ce43b3ae6bd1a27a7233156e32938ee76bd6015a27910a7fa/data'
        : 'community.wave.seqera.io/library/python_pyyaml_faiss-cpu_pip_scib-metrics:f0a647b4acd07c42'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    val max_cells
    val subsample_strategy
    val subsample_seed
    val metric_profile
    val neighbor_backend

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
