process ADATA_READRDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/90/90f55177e533b8e70e61453553c784658ce4a092d6837d18f3aa5ab329ef2e51/data'
        : 'community.wave.seqera.io/library/adata_readrds:13e30e5d7d25206e'}"

    input:
    tuple val(meta), path(rds)
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr'])) {
        error("anndata_format must be 'h5ad' or 'zarr', got '${anndata_format}'")
    }
    template('readrds.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr'])) {
        error("anndata_format must be 'h5ad' or 'zarr', got '${anndata_format}'")
    }
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : "touch ${prefix}.h5ad"
    """
    ${anndata_stub}
    touch versions.yml
    """
}
