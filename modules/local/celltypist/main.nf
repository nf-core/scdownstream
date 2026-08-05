process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39223737b9353409f4db6c18e7143b106ad9f322cd3a82d7741d89e49a5ba2fa/data'
:         'community.wave.seqera.io/library/celltypist:a4fa2621f1596f41' }"

    input:
    tuple val(meta), path(h5ad), val(symbol_col)
    val models

    output:
    tuple val(meta), path("*.h5ad")                  , emit: h5ad
    tuple val(meta), path("*.parquet")                   , emit: obs
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
    touch ${prefix}.parquet
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
