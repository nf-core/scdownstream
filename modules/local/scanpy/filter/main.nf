process SCANPY_FILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4046c811cf4df04aa2f9e7039f452117ea355ba87b2a7505b2c444e3c577e3e/data' :
        'community.wave.seqera.io/library/python_pyyaml_anndata_scanpy:f7441a003c01f440' }"

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
    path "${prefix}_qc_histograms.png"     , emit: plots, optional: true
    path "${prefix}_qc_histograms_mqc.json", emit: multiqc_files, optional: true
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
        touch ${prefix}_qc_histograms.png
        touch ${prefix}_qc_histograms_mqc.json
    fi
    """
}
