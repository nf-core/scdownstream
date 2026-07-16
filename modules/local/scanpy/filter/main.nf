process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/45/45339bf761a2cf0cdb058492bc37f3df8b05b363731d491d1d3a14e9ba0b8f55/data'
        : 'community.wave.seqera.io/library/harmonypy_anndata_leidenalg_numpy_pruned:43066d5f86f18261'}"

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
    path "*.png"                           , emit: plots, optional: true
    path "*_mqc.json"                      , emit: multiqc_files, optional: true
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    section_name = task.ext.section_name ?: "Filter threshold histograms"
    description = task.ext.description ?: "QC metric histograms with applied filter thresholds"
    if ("${prefix}.h5ad" == "${h5ad}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('filter.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    section_name = task.ext.section_name ?: "Filter threshold histograms"
    description = task.ext.description ?: "QC metric histograms with applied filter thresholds"
    """
    touch ${prefix}.h5ad
    touch versions.yml

    if [ "${plot}" = "true" ]; then
        for metric in log1p_total_counts log1p_n_genes_by_counts pct_counts_in_top_20_genes pct_counts_mt pct_counts_ribo pct_counts_hb total_counts n_genes_by_counts; do
            touch ${prefix}_\${metric}.png
            touch ${prefix}_\${metric}_mqc.json
        done
    fi
    """
}
