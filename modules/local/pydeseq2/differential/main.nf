process PYDESEQ2_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b8/b8787e6a23696ee360a647491dcbf90375d67ab0d273f671696636299738ece6/data'
        : 'community.wave.seqera.io/library/pydeseq2_differential:15a563a547b5b3bd' }"

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
