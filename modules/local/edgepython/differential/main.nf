process EDGEPYTHON_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/dd2bdd0dd3afef572c1ae41ca13cc1cd7a3e5ca62410c39934b10964ee916f71/data' :
        'community.wave.seqera.io/library/edgepython_differential:663a87cf9609bfbd' }"

    input:
    tuple val(meta), path(h5ad)
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
