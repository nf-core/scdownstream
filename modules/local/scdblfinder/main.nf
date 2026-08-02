process SCDBLFINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/11/11ce5488f05e7cfea5577fff9bbe172dbf1e3d427ded8fa09f3350a84f4e24ed/data' :
        'community.wave.seqera.io/library/scdblfinder:bf67b6150784b907' }"

    input:
    tuple val(meta), path(h5ad), val(dbr), val(batch_col)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.csv"), emit: predictions
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('scdblfinder.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.csv
    touch versions.yml
    """
}
