/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPROCESS             } from '../subworkflows/local/preprocess'
include { COMBINE                } from '../subworkflows/local/combine'
include { CELLTYPE_ASSIGNMENT    } from '../subworkflows/local/celltype_assignment'
include { CLUSTER                } from '../subworkflows/local/cluster'
include { FINALIZE               } from '../subworkflows/local/finalize'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
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
    ch_base // value channel: [ val(meta), h5ad, scvi_model, model_type ]

    main:

    ch_versions = Channel.empty()
    ch_obs = Channel.empty()
    ch_obsm = Channel.empty()
    ch_obsp = Channel.empty()
    ch_uns = Channel.empty()
    ch_layers = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Per-sample preprocessing
    //

    PREPROCESS(ch_samplesheet)
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS.out.multiqc_files)

    //
    // Combine samples and perform integration
    //

    COMBINE(PREPROCESS.out.h5ad, ch_base)
    ch_versions = ch_versions.mix(COMBINE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COMBINE.out.multiqc_files)
    ch_obs = ch_obs.mix(COMBINE.out.obs)
    ch_obsm = ch_obsm.mix(COMBINE.out.obsm)
    ch_layers = ch_layers.mix(COMBINE.out.layers)

    //
    // Perform automated celltype assignment
    //
    CELLTYPE_ASSIGNMENT(COMBINE.out.h5ad)
    ch_versions = ch_versions.mix(CELLTYPE_ASSIGNMENT.out.versions)
    ch_obs = ch_obs.mix(CELLTYPE_ASSIGNMENT.out.obs)

    //
    // Perform clustering and per-cluster analysis
    //

    CLUSTER(COMBINE.out.integrations)
    ch_versions = ch_versions.mix(CLUSTER.out.versions)
    ch_obs = ch_obs.mix(CLUSTER.out.obs)
    ch_obsm = ch_obsm.mix(CLUSTER.out.obsm)
    ch_obsp = ch_obsp.mix(CLUSTER.out.obsp)
    ch_uns = ch_uns.mix(CLUSTER.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(CLUSTER.out.multiqc_files)

    FINALIZE(COMBINE.out.h5ad, ch_obs, ch_obsm, ch_obsp, ch_uns, ch_layers)
    ch_versions = ch_versions.mix(FINALIZE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
