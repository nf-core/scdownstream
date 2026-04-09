#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/scdownstream
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/scdownstream
    Website: https://nf-co.re/scdownstream
    Slack  : https://nfcore.slack.com/channels/scdownstream
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCDOWNSTREAM            } from './workflows/scdownstream'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scdownstream_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scdownstream_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_SCDOWNSTREAM {

    take:
    samplesheet                   // channel: samplesheet read in from --input
    ch_base                       // channel: [ val(meta), path(h5ad) ]
    is_extension                  //   value: boolean
    ch_input                      //    file: samplesheet.csv
    ambient_correction            //   value: string
    ambient_corrected_integration //   value: boolean
    doublet_detection             //   value: string
    doublet_detection_threshold   //   value: integer
    scvi_max_epochs               //   value: integer
    mito_genes                    //   value: string
    sample_n                      //   value: string
    sample_fraction               //   value: string
    cell_cycle_scoring            //   value: boolean
    s_genes                       //    path: file or []
    g2m_genes                     //    path: file or []
    qc_only                       //   value: boolean
    celldex_reference             //   value: string
    celltypist_model              //   value: string
    unify_gene_symbols            //   value: boolean
    duplicate_var_resolution      //   value: string
    aggregate_isoforms            //   value: boolean
    integration_hvgs              //   value: integer
    integration_methods           //   value: string
    integration_excluded_genes    //   value: string
    scvi_model                    //   value: string
    scanvi_model                  //   value: string
    scvi_categorical_covariates   //   value: string
    scvi_continuous_covariates    //   value: string
    scimilarity_model             //   value: string
    skip_liana                    //   value: boolean
    skip_rankgenesgroups          //   value: boolean
    base_embeddings               //   value: string
    base_condition_col            //   value: string
    base_label_col                //   value: string
    cluster_per_label             //   value: boolean
    cluster_global                //   value: boolean
    clustering_resolutions        //   value: string
    pseudobulk                    //   value: boolean
    pseudobulk_groupby_labels     //   value: string
    pseudobulk_min_num_cells      //   value: integer
    prep_cellxgene                //   value: boolean
    outdir                        //   value: string
    multiqc_config                //   value: string
    multiqc_logo                  //   value: string
    multiqc_methods_description   //   value: string

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SCDOWNSTREAM (
        samplesheet,
        ch_base,
        is_extension,
        ch_input,
        ambient_correction,
        ambient_corrected_integration,
        doublet_detection,
        doublet_detection_threshold,
        scvi_max_epochs,
        mito_genes,
        sample_n,
        sample_fraction,
        cell_cycle_scoring,
        s_genes,
        g2m_genes,
        qc_only,
        celldex_reference,
        celltypist_model,
        unify_gene_symbols,
        duplicate_var_resolution,
        aggregate_isoforms,
        integration_hvgs,
        integration_methods,
        integration_excluded_genes,
        scvi_model,
        scanvi_model,
        scvi_categorical_covariates,
        scvi_continuous_covariates,
        scimilarity_model,
        skip_liana,
        skip_rankgenesgroups,
        base_embeddings,
        base_label_col,
        base_condition_col,
        cluster_per_label,
        cluster_global,
        clustering_resolutions,
        pseudobulk,
        pseudobulk_groupby_labels,
        pseudobulk_min_num_cells,
        prep_cellxgene,
        outdir,
        multiqc_config,
        multiqc_logo,
        multiqc_methods_description
    )
    emit:
    multiqc_report = SCDOWNSTREAM.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    ch_base_adata = params.base_adata
            ? channel.value([[id: "base"], file(params.base_adata, checkIfExists: true)])
            : channel.value([[], []])

    def s_genes_file = params.cell_cycle_scoring
        ? file(params.s_genes ?: "${projectDir}/assets/cell_cycle_genes/${params.species}_s_genes.txt", checkIfExists: true)
        : []
    def g2m_genes_file = params.cell_cycle_scoring
        ? file(params.g2m_genes ?: "${projectDir}/assets/cell_cycle_genes/${params.species}_g2m_genes.txt", checkIfExists: true)
        : []

    NFCORE_SCDOWNSTREAM (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_base_adata,
        params.base_adata != null,
        params.input,
        params.ambient_correction,
        params.ambient_corrected_integration,
        params.doublet_detection,
        params.doublet_detection_threshold,
        params.scvi_max_epochs,
        params.mito_genes,
        params.sample_n,
        params.sample_fraction,
        params.cell_cycle_scoring,
        s_genes_file,
        g2m_genes_file,
        params.qc_only,
        params.celldex_reference,
        params.celltypist_model,
        params.unify_gene_symbols,
        params.duplicate_var_resolution,
        params.aggregate_isoforms,
        params.integration_hvgs,
        params.integration_methods,
        params.integration_excluded_genes,
        params.scvi_model,
        params.scanvi_model,
        params.scvi_categorical_covariates,
        params.scvi_continuous_covariates,
        params.scimilarity_model,
        params.skip_liana,
        params.skip_rankgenesgroups,
        params.base_embeddings,
        params.base_label_col,
        params.base_condition_col,
        params.cluster_per_label,
        params.cluster_global,
        params.clustering_resolutions,
        params.pseudobulk,
        params.pseudobulk_groupby_labels,
        params.pseudobulk_min_num_cells,
        params.prep_cellxgene,
        params.outdir,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_SCDOWNSTREAM.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
