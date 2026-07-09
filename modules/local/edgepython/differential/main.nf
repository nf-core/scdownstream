process EDGEPYTHON_DIFFERENTIAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ef/eff49ff9ce5ae24018eb162cdd5b32761901b28514fff1f61781ac1825a44c61/data'
        : 'community.wave.seqera.io/library/edgepython_pseudobulk:d8a69e7fb7ac4371' }"

    input:
    tuple val(meta), path(h5ad)
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
