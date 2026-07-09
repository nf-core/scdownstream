//
// Subworkflow with functionality specific to the nf-core/scdownstream pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/scdownstream ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/scdownstream/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    ch_samplesheet = params.input
        ? channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")).map {
            sample -> validateInputSamplesheet(sample)
        }
        : channel.empty()

    emit:
    samplesheet        = ch_samplesheet
    versions           = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Default analysis plan: one row with empty wildcards (match all clusterings, all analyses)
//
def analysisPlanToList() {
    params.analysis_plan
        ? samplesheetToList(params.analysis_plan, "${projectDir}/assets/schema_analysis_plan.json")
            .collect { row -> row[0] }
        : [[integration: null, subset: null, resolution: null, analyses: null, de_methods: null]]
}

def matchesAnalysisPlanRow(row, meta, resolution = null) {
    (!row.integration || row.integration == meta.integration) &&
    (!row.subset || row.subset == meta.subset) &&
    (resolution == null || !row.resolution || (row.resolution as String) == (resolution as String))
}

def matchingAnalysisPlanRows(rows, meta, resolution = null) {
    rows.findAll { row -> matchesAnalysisPlanRow(row, meta, resolution) }
}

def analysesFromPlanRows(rows) {
    if (!rows || rows.any { row -> !row.analyses }) {
        return [:]
    }
    [
        analyses: rows
            .collectMany { row -> row.analyses.split(',').collect { token -> token.trim() } }
            .toSet(),
    ]
}

def deMethodsFromPlanRows(rows) {
    def methods = rows
        .findAll { row -> row.de_methods }
        .collectMany { row -> row.de_methods.split(',').collect { token -> token.trim() } }
        .findAll { method -> method }
        .toSet()
    if (methods.isEmpty()) {
        return [:]
    }
    [de_methods: methods]
}

def resolveDeMethods(meta, default_methods) {
    if (meta.de_methods) {
        return meta.de_methods as Set
    }
    return default_methods.split(',').collect { token -> token.trim() }.findAll { token -> token } as Set
}

def cytetypeEligible(meta, cytetype_study_context) {
    cytetype_study_context && (meta.analyses == null || 'cytetype' in meta.analyses)
}

def resolveDeMethodsWithPrerequisites(meta, default_methods, cytetype_study_context = '') {
    def resolved = resolveDeMethods(meta, default_methods)
    def extra = [:]
    if (cytetypeEligible(meta, cytetype_study_context)) {
        if (!('wilcoxon' in resolved)) {
            resolved = resolved + 'wilcoxon'
        }
        extra.cytetype_prerequisite = true
    }
    return [de_methods_resolved: resolved] + extra
}

def rankGenesGroupsAnalysisEnabled(meta) {
    meta.analyses == null || 'de' in meta.analyses || meta.cytetype_prerequisite
}

def rankGenesGroupsMethods() {
    ['wilcoxon', 't-test', 't-test_overestim_var', 'logreg'] as Set
}

def rankGenesGroupsMethodSlug(method) {
    method.replace('-', '_')
}

def rankGenesGroupsMethodKey(method, single_method = false) {
    single_method ? 'rank_genes_groups' : "rank_genes_groups_${rankGenesGroupsMethodSlug(method)}"
}

def preferredRankGenesGroupsMethod(methods) {
    def order = ['wilcoxon', 't-test', 't-test_overestim_var', 'logreg']
    def resolved = methods.intersect(rankGenesGroupsMethods())
    order.find { method -> method in resolved } ?: resolved.toList()[0]
}

def validDeMethods() {
    (rankGenesGroupsMethods() + ['pydeseq2', 'edgepython', 'edgepython_sc']) as Set
}

def cellLevelDeMethods() {
    (rankGenesGroupsMethods() + ['edgepython_sc']) as Set
}

def pseudobulkDeMethods() {
    ['pydeseq2', 'edgepython'] as Set
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    if (!params.input && !(params.base_adata && params.base_label_col && (params.base_embeddings || params.integrate_per_label))) {
        throw new Exception("Either an input samplesheet or (base_adata && base_label_col && (base_embeddings || integrate_per_label)) must be provided")
    }

    if (params.qc_only && !params.input) {
        throw new Exception("If qc_only is set to true, an input samplesheet must be provided")
    }

    if (params.integrate_per_label_whitelist && !params.integrate_per_label) {
        throw new Exception("integrate_per_label_whitelist requires integrate_per_label to be true")
    }

    def integration_methods = params.integration_methods.split(',').collect { it -> it.trim().toLowerCase() }
    def is_extension = params.input && params.base_adata
    def is_per_label_base_integration = !params.input && params.base_adata && params.integrate_per_label

    if (is_extension && (integration_methods - ['scvi', 'scanvi', 'scimilarity', 'symphony']).size() > 0) {
        throw new Exception("Only scvi, scanvi, scimilarity and symphony integration methods are supported if base_adata is provided")
    }

    if (is_extension && 'scvi' in integration_methods && !params.scvi_model) {
        throw new Exception("If base_adata is provided and scvi is used as integration method, scvi_model must be provided.")
    }

    if (is_extension && 'scanvi' in integration_methods && !params.scanvi_model) {
        throw new Exception("If base_adata is provided and scanvi is used as integration method, scanvi_model must be provided.")
    }

    if ((is_extension || is_per_label_base_integration) && 'scimilarity' in integration_methods && !params.scimilarity_model) {
        throw new Exception("If base_adata is provided and scimilarity is used as integration method, scimilarity_model must be provided.")
    }

    if (is_extension && 'symphony' in integration_methods && !params.symphony_reference) {
        throw new Exception("If base_adata is provided and symphony is used as integration method, symphony_reference must be provided.")
    }

    // Validate sample_n and sample_fraction parameters
    if (params.sample_n && params.sample_fraction) {
        throw new Exception("Both sample_n and sample_fraction are set. Please use only one of them.")
    }

    def de_methods = params.de_methods.split(',').collect { token -> token.trim() }.findAll { token -> token }
    def invalid_de_methods = de_methods.findAll { method -> !(method in validDeMethods()) }
    if (invalid_de_methods) {
        throw new Exception("Invalid de_methods: ${invalid_de_methods.join(', ')}. Valid options: ${validDeMethods().join(', ')}")
    }

    def pseudobulk_methods = de_methods.intersect(pseudobulkDeMethods() as List)
    if (pseudobulk_methods) {
        def plan_rows = params.analysis_plan ? analysisPlanToList() : []
        def plan_has_pseudobulk_de = plan_rows.any { row -> row.analyses && 'pseudobulk_de' in row.analyses.split(',').collect { token -> token.trim() } }
        if (!plan_has_pseudobulk_de) {
            log.warn "de_methods includes ${pseudobulk_methods.join(', ')} but no analysis_plan row lists pseudobulk_de in analyses. Pseudobulk DE will not run."
        }
    }

    def plan_rows_with_pb_methods = params.analysis_plan
        ? analysisPlanToList().findAll { row ->
            row.de_methods && (row.de_methods.split(',').collect { token -> token.trim() }.intersect(pseudobulkDeMethods() as List))
        }
        : []
    if (plan_rows_with_pb_methods) {
        def missing_token = plan_rows_with_pb_methods.findAll { row ->
            !row.analyses || !('pseudobulk_de' in row.analyses.split(',').collect { token -> token.trim() })
        }
        if (missing_token) {
            log.warn "analysis_plan rows specify pseudobulk de_methods (${pseudobulkDeMethods().join(', ')}) without pseudobulk_de in analyses. Those engines will not run for matching clusterings."
        }
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (meta, filtered, unfiltered) = input
    if (!filtered && !unfiltered) {
        throw new Exception("Both filtered and unfiltered files are missing for sample ${meta.id}")
    }

    return input
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
