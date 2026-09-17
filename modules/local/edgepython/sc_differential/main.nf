process EDGEPYTHON_SCDIFFERENTIAL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/387ee8bae78078479f4247c84d00dd1a6f97ba49315f807d953d409b19be9823/data' :
        'community.wave.seqera.io/library/edgepython_sc_differential:2fe5a778fe82c638' }"

    input:
    tuple val(meta), path(h5ad)
    val(donor_col)
    val(condition_col)
    val(celltype_col)
    val(celltype_value)
    val(reference_condition)
    path interesting_genes

    output:
    tuple val(meta), path("${prefix}_*_results.csv"), emit: results
    path "*.png"                                    , emit: plots, optional: true
    path "*_mqc.json"                               , emit: multiqc_files, optional: true
    path "versions.yml"                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('sc_differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_treatment_results.csv
    touch ${prefix}_treatment_volcano.png
    touch ${prefix}_treatment_volcano_mqc.json
    touch versions.yml
    """
}
