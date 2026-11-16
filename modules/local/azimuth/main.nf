process CELLTYPES_AZIMUTH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fefe5f28dd65978dde0b537e9d7f4c83a1fa173e73c5fcd4ce53badc01f64a5b/data'
        : 'community.wave.seqera.io/library/anndata_numpy_pandas_pyarrow_pruned:09e83165233226e4'}"

    input:
    tuple val(meta), path(h5ad), val(symbol_col), val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}.parquet"), emit: obs
    tuple val(meta), path("X_azimuth.parquet"), emit: obsm
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    counts_layer = counts_layer ?: "X"
    template('azimuth.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch X_azimuth.parquet
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
