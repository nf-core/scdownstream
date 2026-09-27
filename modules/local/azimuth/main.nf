process CELLTYPES_AZIMUTH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/3274eab1b78962b7144695e36e24a902dd1117d3f266a402ac3edbce4b51836d/data'
:         'community.wave.seqera.io/library/anndata_numpy_pandas_python_pruned:7e74e77ee7540c2e' }"

    input:
    tuple val(meta), path(h5ad), val(symbol_col), val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad")          , emit: h5ad
    tuple val(meta), path("${prefix}.pkl")           , emit: obs
    tuple val(meta), path("X_azimuth.pkl")           , emit: obsm
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    path "versions.yml"                              , emit: versions, topic: versions

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
    touch ${prefix}.pkl
    touch X_azimuth.pkl
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
