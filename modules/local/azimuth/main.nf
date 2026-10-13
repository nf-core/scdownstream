process CELLTYPES_AZIMUTH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1d/1dfc349185b35fbb30c64f1a8dfebd4773e4bbc54f01be2d3006172bab3cd70f/data'
        : 'community.wave.seqera.io/library/anndata_numpy_pandas_python_pruned:82d6d447e5f8fd2b'}"

    input:
    tuple val(meta), path(anndata), val(symbol_col), val(counts_layer)
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    tuple val(meta), path("${prefix}.pkl"), emit: obs
    tuple val(meta), path("X_azimuth.pkl"), emit: obsm
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    counts_layer = counts_layer ?: "X"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (anndata_format && "${prefix}.${anndata_format}" == "${anndata}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('azimuth.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch ${prefix}.pkl
    touch X_azimuth.pkl
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
