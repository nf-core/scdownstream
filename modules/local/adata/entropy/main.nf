process ADATA_ENTROPY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04307b777a83d13f004023a046a677adc1c444edcae99a1232cc5b82e88968df/data'
        : 'community.wave.seqera.io/library/leidenalg_pyarrow_python-igraph_python_pruned:085892d5c5a30e57'}"

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
