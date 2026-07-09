process EDGEPYTHON_SCDIFFERENTIAL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7750a9d7aa07f2f73875f1f81ea2351480496f24bc782ca168364fd3f6747e5/data'
        : 'community.wave.seqera.io/library/edgepython_sc:76c19aa233370b64' }"

    input:
    tuple val(meta), path(h5ad)
    val(donor_col)
    val(condition_col)
    val(celltype_col)
    val(celltype_value)
    val(reference_condition)

    output:
    tuple val(meta), path("${prefix}_results.csv"), emit: results
    path "versions.yml"                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('sc_differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_results.csv
    touch versions.yml
    """
}
