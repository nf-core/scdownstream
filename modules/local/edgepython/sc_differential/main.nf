process EDGEPYTHON_SCDIFFERENTIAL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7967cbf85b1f33815b4e84540f1ba3486dc6e7a87283b5ea5e5fe10391466851/data':
        'community.wave.seqera.io/library/edgepython_sc_differential:d3364fef9746874f' }"

    input:
    tuple val(meta), path(h5ad)
    val(donor_col)
    val(condition_col)
    val(celltype_col)
    val(celltype_value)
    val(reference_condition)

    output:
    tuple val(meta), path("${prefix}_*_results.parquet"), emit: results
    path "versions.yml"                                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('sc_differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_treatment_results.parquet
    touch versions.yml
    """
}
