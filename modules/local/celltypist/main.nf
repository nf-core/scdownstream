process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/88d5c1e86b7599bee9d80818628a4d9dc5de36b9e0915522b170919bc742bcd2/data' :
        'community.wave.seqera.io/library/celltypist:34649e0085cf7125' }"

    input:
    tuple val(meta), path(h5ad), val(symbol_col)
    val models

    output:
    tuple val(meta), path("*.h5ad")                  , emit: h5ad
    tuple val(meta), path("*.pkl")                   , emit: obs
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('celltypist.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.pkl
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
