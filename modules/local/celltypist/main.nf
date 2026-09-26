process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb0b40342ced5d0b963803016966c024757ffa84bb3514f20347bbc1de113c97/data'
:         'community.wave.seqera.io/library/celltypist:5c8a1d7b51207f8e' }"

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
