process EDGEPYTHON_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bc5e56de4a965746ba8d70e4d00e73f71815b799ac6747b0e5608cbbc8ca7623/data':
        'community.wave.seqera.io/library/edgepython_differential:2240484d16e01bd0' }"

    input:
    tuple val(meta), path(h5ad)
    val(reference_condition)

    output:
    tuple val(meta), path("${prefix}_*_results.parquet"), emit: results
    path "versions.yml"                                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stratum_results.parquet
    touch versions.yml
    """
}
