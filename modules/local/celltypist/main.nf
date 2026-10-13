process CELLTYPES_CELLTYPIST {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/374282d1a3622b00af7c63491aba54a7c0a51fba3421efa9b8281734ba67cbd8/data'
        : 'community.wave.seqera.io/library/celltypist:a91a353b6f0c32c6'}"

    input:
    tuple val(meta), path(anndata), val(symbol_col)
    val models
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    tuple val(meta), path("*.pkl"), emit: obs
    tuple val(meta), path("*_annotation_columns.csv"), emit: annotation_columns
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (anndata_format && "${prefix}.${anndata_format}" == "${anndata}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('celltypist.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch ${prefix}.pkl
    echo "obs_column,aggregatable" > ${prefix}_annotation_columns.csv
    touch versions.yml
    """
}
