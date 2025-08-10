/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LOAD_H5AD                            } from '../subworkflows/local/load_h5ad'
include { QUALITY_CONTROL                      } from '../subworkflows/local/quality_control'
include { CELLTYPE_ASSIGNMENT                  } from '../subworkflows/local/celltype_assignment'
include { ADATA_EXTEND as FINALIZE_QC_ANNDATAS } from '../modules/local/adata/extend'
include { COMBINE                              } from '../subworkflows/local/combine'
include { ADATA_SPLITEMBEDDINGS                } from '../modules/local/adata/splitembeddings'
include { CLUSTER                              } from '../subworkflows/local/cluster'
include { PSEUDOBULKING                        } from '../subworkflows/local/pseudobulking'
include { PER_GROUP                            } from '../subworkflows/local/per_group'
include { FINALIZE                             } from '../subworkflows/local/finalize'
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText               } from '../subworkflows/local/utils_nfcore_scdownstream_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCDOWNSTREAM {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_base        // channel: [ val(meta), path(h5ad) ]

    main:

    ch_versions = channel.empty()
    ch_integrations = channel.empty()
    ch_obs = channel.empty()
    ch_var = channel.empty()
    ch_obsm = channel.empty()
    ch_obsp = channel.empty()
    ch_uns = channel.empty()
    ch_layers = channel.empty()
    ch_multiqc_files = channel.empty()

    if (params.input) {
        ch_obs_per_sample = channel.empty()
        ch_var_per_sample = channel.empty()
        ch_obsm_per_sample = channel.empty()
        ch_obsp_per_sample = channel.empty()
        ch_uns_per_sample = channel.empty()
        ch_layers_per_sample = channel.empty()

        //
        // Load/Convert input to h5ad
        //
        LOAD_H5AD(ch_samplesheet)
        ch_h5ad = LOAD_H5AD.out.h5ad
        ch_versions = ch_versions.mix(LOAD_H5AD.out.versions)

        //
        // Quality control per sample
        //
        QUALITY_CONTROL(
            ch_h5ad,
            params.ambient_correction,
            (!params.doublet_detection || params.doublet_detection == 'none') ? [] : params.doublet_detection.split(',').collect { it -> it.trim().toLowerCase() },
            params.doublet_detection_threshold,
            params.mito_genes,
        )
        ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.multiqc_files)
        ch_h5ad = QUALITY_CONTROL.out.h5ad

        //
        // Perform automated celltype assignment
        //
        CELLTYPE_ASSIGNMENT(ch_h5ad.map { meta, h5ad -> [meta, h5ad, meta.symbol_col] })
        ch_versions = ch_versions.mix(CELLTYPE_ASSIGNMENT.out.versions)
        ch_obs_per_sample = ch_obs_per_sample.mix(CELLTYPE_ASSIGNMENT.out.obs)

        FINALIZE_QC_ANNDATAS(
            ch_h5ad.join(ch_obs_per_sample.groupTuple(), remainder: true).join(ch_var_per_sample.groupTuple(), remainder: true).join(ch_obsm_per_sample.groupTuple(), remainder: true).join(ch_obsp_per_sample.groupTuple(), remainder: true).join(ch_uns_per_sample.groupTuple(), remainder: true).join(ch_layers_per_sample.groupTuple(), remainder: true).map { meta, h5ad, obs, var, obsm, obsp, uns, layers ->
                [meta, h5ad, obs ?: [], var ?: [], obsm ?: [], obsp ?: [], uns ?: [], layers ?: []]
            }
        )
        ch_h5ad = FINALIZE_QC_ANNDATAS.out.h5ad
        ch_versions = ch_versions.mix(FINALIZE_QC_ANNDATAS.out.versions)

        if (!params.qc_only) {
            //
            // Combine samples and perform integration
            //
            COMBINE(ch_h5ad, ch_base)
            ch_versions = ch_versions.mix(COMBINE.out.versions)
            ch_obs = ch_obs.mix(COMBINE.out.obs)
            ch_var = ch_var.mix(COMBINE.out.var)
            ch_obsm = ch_obsm.mix(COMBINE.out.obsm)
            ch_integrations = ch_integrations.mix(COMBINE.out.integrations)
            ch_finalization_base = COMBINE.out.h5ad

            ch_label_grouping = COMBINE.out.h5ad_inner
            grouping_col = "label"
        }
    }
    else {
        ch_embeddings = channel.value(params.base_embeddings.split(',').collect { it -> it.trim() })

        ADATA_SPLITEMBEDDINGS(ch_base, ch_embeddings)
        ch_versions = ch_versions.mix(ADATA_SPLITEMBEDDINGS.out.versions)
        ch_integrations = ch_integrations.mix(
            ADATA_SPLITEMBEDDINGS.out.h5ad.map { _meta, h5ads -> h5ads }.flatten().map { h5ad -> [[id: h5ad.simpleName, integration: h5ad.simpleName], h5ad] }
        )

        ch_finalization_base = ch_base
        ch_label_grouping = ch_base
        grouping_col = params.base_label_col
    }

    //
    // Perform clustering and per-cluster analysis
    //
    if (!params.qc_only) {
        CLUSTER(
            ch_integrations,
            params.cluster_per_label,
            params.cluster_global,
            params.input ? "label" : params.base_label_col,
            params.clustering_resolutions.split(','),
            "batch",
            "X_emb",
        )
        ch_versions = ch_versions.mix(CLUSTER.out.versions)
        ch_obs = ch_obs.mix(CLUSTER.out.obs)
        ch_obsm = ch_obsm.mix(CLUSTER.out.obsm)
        ch_multiqc_files = ch_multiqc_files.mix(CLUSTER.out.multiqc_files)

        if (params.pseudobulk) {
            PSEUDOBULKING(
                CLUSTER.out.h5ad_clustering,
                params.pseudobulk_groupby_labels.split(','),
                params.pseudobulk_min_num_cells,
                "X",
            )
            ch_versions = ch_versions.mix(PSEUDOBULKING.out.versions)
        }

        ch_h5ad_both = CLUSTER.out.h5ad_clustering.map { meta, h5ad -> [meta + [obs_key: "${meta.id}_leiden"], h5ad] }

        PER_GROUP(
            // Run on each clustering resolution for each embedding
            ch_h5ad_both.mix(
                // And on the label column for each embedding
                CLUSTER.out.h5ad_neighbors.map { meta, h5ad -> [meta + [obs_key: grouping_col], h5ad] }
            ),
            // Run on each clustering (there is one clustering per embedding and resolution)
            ch_h5ad_both.mix(
                // And on the label column
                ch_label_grouping.map { meta, h5ad -> [meta + [obs_key: grouping_col], h5ad] }
            ),
        )
        ch_versions = ch_versions.mix(PER_GROUP.out.versions)
        ch_uns = ch_uns.mix(PER_GROUP.out.uns)
        ch_multiqc_files = ch_multiqc_files.mix(PER_GROUP.out.multiqc_files)

        FINALIZE(ch_finalization_base, ch_obs, ch_var, ch_obsm, ch_obsp, ch_uns, ch_layers)
        ch_versions = ch_versions.mix(FINALIZE.out.versions)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'scdownstream_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
