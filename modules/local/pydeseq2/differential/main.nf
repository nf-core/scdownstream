process PYDESEQ2_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/07aa2f42d62525549d9a7265a8678c794e87b3ac403514fa1065b879e9e6cffb/data'
        : 'community.wave.seqera.io/library/pydeseq2_differential:4c8249b844ab7a1a' }"

    input:
    tuple val(meta), path(h5ad)
    val(design_formula)
    val(reference_condition)

    output:
    tuple val(meta), path("${prefix}_*_results.csv"), emit: results
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stratum_results.csv
    touch versions.yml
    """
}
