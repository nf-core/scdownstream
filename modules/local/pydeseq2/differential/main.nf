process PYDESEQ2_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/779736ea3b01e9f685d5a23697524e4a558e59284a73a51c779713cb697a4f59/data':
        'community.wave.seqera.io/library/pydeseq2_differential:42831fea4fd7c416' }"

    input:
    tuple val(meta), path(h5ad)
    val(design_formula)
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
