process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c43892b5ea2991c98bd69e6da333b902a1e54635e24ca2cf3fd4c5972efb4384/data'
        : 'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:e5a54bd2b6c9720d'}"

    input:
    tuple val(meta), path(h5ad)
    val symbol_col
    val min_genes
    val min_cells
    val min_counts_gene
    val min_counts_cell
    val max_mito_percentage
    val min_ribo_percentage
    val max_hb_percentage
    val log1p_total_counts_nmads
    val log1p_n_genes_by_counts_nmads
    val pct_counts_in_top_20_genes_nmads
    val pct_counts_mt_nmads
    path mito_genes
    val plot

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "${prefix}_qc_histograms.png", emit: plots, optional: true
    path "${prefix}_qc_histograms_mqc.json", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('filter.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml

    if [ "${plot}" = "true" ]; then
        touch ${prefix}_qc_histograms.png
        touch ${prefix}_qc_histograms_mqc.json
    fi
    """
}
