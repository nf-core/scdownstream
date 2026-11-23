process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0c5ec8266eda92dde1726fbd2bcecde46e376eee341f0388c0134fcdcfe92ce/data'
        : 'community.wave.seqera.io/library/celltypist:ad6d33f24b6d33ff'}"

    input:
    tuple val(meta), path(h5ad), val(symbol_col)
    val models
    path model_files, arity: '0..*', stageAs: 'celltypist_models/*'

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.parquet"), emit: obs
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    path "versions.yml", emit: versions, topic: versions

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
