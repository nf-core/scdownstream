process PYDESEQ2_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/ded0dc26fe39d71f6eb1d888077b1abfa9fca33bd36fb51c295c5429abb4dc04/data' :
        'community.wave.seqera.io/library/pydeseq2_differential:c83f772671389253' }"

    input:
    tuple val(meta), path(h5ad)
    val(design_formula)
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
    template('differential.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stratum_results.csv
    touch ${prefix}_stratum_volcano.png
    touch ${prefix}_stratum_volcano_mqc.json
    touch versions.yml
    """
}
