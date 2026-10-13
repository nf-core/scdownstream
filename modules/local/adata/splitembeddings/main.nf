process ADATA_SPLITEMBEDDINGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf'}"

    input:
    tuple val(meta), path(anndata)
    val embeddings
    val anndata_format

    output:
    tuple val(meta), path("*.{h5ad,zarr}"), emit: anndata
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (!anndata_format) {
        error("ADATA_SPLITEMBEDDINGS only outputs the AnnData object; anndata_format must not be empty")
    }
    template('split_embeddings.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def stub_cmd = anndata_format == 'zarr'
        ? "mkdir scvi.zarr && touch scvi.zarr/zarr.json"
        : "touch scvi.h5ad"
    """
    ${stub_cmd}
    touch versions.yml
    """
}
