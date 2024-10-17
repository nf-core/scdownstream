/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS             } from '../subworkflows/local/preprocess'
include { CELLTYPE_ASSIGNMENT    } from '../subworkflows/local/celltype_assignment'
include { COMBINE                } from '../subworkflows/local/combine'
include { ADATA_SPLITEMBEDDINGS  } from '../modules/local/adata/splitembeddings'
include { CLUSTER                } from '../subworkflows/local/cluster'
include { FINALIZE               } from '../subworkflows/local/finalize'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scdownstream_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCDOWNSTREAM {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_base  // value channel: [ val(meta), path(h5ad) ]
    ch_reference_model // value channel: [ val(meta), path(model) ]

    main:

    ch_versions = Channel.empty()
    ch_integrations = Channel.empty()
    ch_obs = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obsp = Channel.empty()
    ch_uns = Channel.empty()
    ch_layers = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (params.input) {
        //
        // Per-sample preprocessing
        //

        PREPROCESS(ch_samplesheet)
        ch_versions = ch_versions.mix(PREPROCESS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS.out.multiqc_files)
        ch_h5ad = PREPROCESS.out.h5ad

        //
        // Perform automated celltype assignment
        //

        if (!params.preprocess_only) {
            CELLTYPE_ASSIGNMENT(ch_h5ad)
            ch_versions = ch_versions.mix(CELLTYPE_ASSIGNMENT.out.versions)
            ch_h5ad = CELLTYPE_ASSIGNMENT.out.h5ad

            //
            // Combine samples and perform integration
            //

            COMBINE(ch_h5ad, ch_base, ch_reference_model)
            ch_versions      = ch_versions.mix(COMBINE.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(COMBINE.out.multiqc_files)
            ch_obs           = ch_obs.mix(COMBINE.out.obs)
            ch_obsm          = ch_obsm.mix(COMBINE.out.obsm)
            ch_layers        = ch_layers.mix(COMBINE.out.layers)
            ch_integrations  = ch_integrations.mix(COMBINE.out.integrations)

            ch_finalization_base = COMBINE.out.h5ad
        }
    } else {
        ch_embeddings = Channel.value(params.base_embeddings.split(',').collect{it.trim()})

        ADATA_SPLITEMBEDDINGS(ch_base, ch_embeddings)
        ch_versions = ch_versions.mix(ADATA_SPLITEMBEDDINGS.out.versions)
        ch_integrations = ch_integrations.mix(
            ADATA_SPLITEMBEDDINGS.out.h5ad
                .map{meta, h5ads -> h5ads}
                .flatten()
                .map{h5ad -> [[id: h5ad.simpleName, integration: h5ad.simpleName], h5ad]}
        )

        ch_finalization_base = ch_base
    }

    //
    // Perform clustering and per-cluster analysis
    //

    if (!params.preprocess_only) {
        CLUSTER(ch_integrations)
        ch_versions = ch_versions.mix(CLUSTER.out.versions)
        ch_obs = ch_obs.mix(CLUSTER.out.obs)
        ch_obsm = ch_obsm.mix(CLUSTER.out.obsm)
        ch_obsp = ch_obsp.mix(CLUSTER.out.obsp)
        ch_uns = ch_uns.mix(CLUSTER.out.uns)
        ch_multiqc_files = ch_multiqc_files.mix(CLUSTER.out.multiqc_files)

        FINALIZE(ch_finalization_base, ch_obs, ch_obsm, ch_obsp, ch_uns, ch_layers)
        ch_versions = ch_versions.mix(FINALIZE.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
