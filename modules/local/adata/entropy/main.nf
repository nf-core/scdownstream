process ADATA_ENTROPY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8f/8f8e43bfea40f384db72a635038d6dfa776155850b820ae063a3ab458957b561/data'
        : 'community.wave.seqera.io/library/leidenalg_pyarrow_python-igraph_python_pruned:8631c4d5e6493952'}"

    input:
    tuple val(meta), path(h5ad)
    val group_col
    val entropy_col
    val plot_basis

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}.parquet", emit: obs
    path "${prefix}.png", emit: plots, optional: true
    path "${prefix}_mqc.json", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    template('entropy.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet

    if [ ${plot_basis ? 'true' : 'false'} ]; then
        touch ${prefix}.png
        touch ${prefix}_mqc.json
    fi

    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
