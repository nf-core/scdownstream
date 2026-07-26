process EDGEPYTHON_SCDIFFERENTIAL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/50/5061ec6a602be93926c3bb5a5e1f760d93a599b6a165a3c902d60de78bfc9953/data'
        : 'community.wave.seqera.io/library/edgepython_sc_differential:ebc200814cba178c' }"

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
