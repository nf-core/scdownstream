process SCANPY_SCRUBLET {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f553feb2c9ac1149defe1f6307fd8921872957e8b0f4b31751920d803da2471/data'
        : 'community.wave.seqera.io/library/pyarrow_pynndescent_python_pyyaml_pruned:c3db4c55d39c6c42'}"

    input:
    tuple val(meta), path(h5ad)
    val batch_col

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.parquet"), emit: predictions
    tuple val(meta), path("${prefix}.scores.parquet"), emit: scores
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_scrublet"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('scrublet.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_scrublet"
    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    """
    # These are needed to prevent errors during import of scanpy
    # when using singularity/apptainer
    export MPLCONFIGDIR=./tmp/mpl
    export NUMBA_CACHE_DIR=./tmp/numba

    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch ${prefix}.scores.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        scanpy: \$(python3 -c 'import scanpy; print(scanpy.__version__)')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
    END_VERSIONS
    """
}
