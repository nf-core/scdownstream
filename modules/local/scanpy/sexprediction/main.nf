process SCANPY_SEXPREDICTION {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f8e43bfea40f384db72a635038d6dfa776155850b820ae063a3ab458957b561/data'
        : 'community.wave.seqera.io/library/leidenalg_pyarrow_python-igraph_python_pruned:8631c4d5e6493952'}"

    input:
    tuple val(meta), path(h5ad)
    path y_genes
    path x_genes
    val symbol_col

    output:
    tuple val(meta), path("${prefix}.parquet"), emit: obs
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('sexprediction.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.parquet
    touch ${prefix}.h5ad
    touch ${prefix}.tsv
    touch versions.yml
    """
}
