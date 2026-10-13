process ADATA_SPLITCOL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/de2455bdcbc53723a8d942ffb0dec7559690ef4bea6b380fea2edba7c5d637b9/data'
        : 'community.wave.seqera.io/library/anndata_python_pyyaml:e77bb4f393d60a08'}"

    input:
    tuple val(meta), path(anndata)
    val column
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
        error("ADATA_SPLITCOL only outputs the AnnData object; anndata_format must not be empty")
    }
    template('split_column.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def stub_cmd = anndata_format == 'zarr'
        ? "mkdir A.zarr B.zarr && touch A.zarr/zarr.json B.zarr/zarr.json"
        : "touch A.h5ad B.h5ad"
    """
    ${stub_cmd}
    touch versions.yml
    """
}
