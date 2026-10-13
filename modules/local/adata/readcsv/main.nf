process ADATA_READCSV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0d347e0d943cb8d0624a8c73247db0c8e96e7c130fe0fb3125aeb7d1b140e79/data'
        : 'community.wave.seqera.io/library/anndata_python_pyyaml:5ac0fd4d280528a5'}"

    input:
    tuple val(meta), path(csv)
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (!anndata_format) {
        error("ADATA_READCSV only outputs the AnnData object; anndata_format must not be empty")
    }
    template('readcsv.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch versions.yml
    """
}
