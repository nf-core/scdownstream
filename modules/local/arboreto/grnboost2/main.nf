process ARBORETO_GRNBOOST2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e04242b1d2b87f2232fa0e7edd211a2c466fac1f6fc4cdcae47d613e7ea35d4/data'
        : 'community.wave.seqera.io/library/arboreto_anndata_pyarrow_python_pruned:b1b1c656e2792807'}"

    input:
    tuple val(meta), path(h5ad)
    path tfs
    val n_top_genes
    val min_cells
    val max_links
    val seed

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    path "${prefix}.parquet", emit: uns
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    obs_key = meta.obs_key ?: "leiden"
    template('grnboost2.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.parquet
    touch versions.yml
    """
}
