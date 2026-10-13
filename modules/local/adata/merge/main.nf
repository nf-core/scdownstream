process ADATA_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf'}"

    input:
    tuple val(meta), path(anndata, stageAs: 'input/?/*')
    tuple val(meta2), path(base)
    val force_obs_cols
    val anndata_format

    output:
    tuple val(meta), path("${prefix}_outer.{h5ad,zarr}"), emit: outer
    tuple val(meta), path("${prefix}_integrate.{h5ad,zarr}"), emit: integrate
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (!anndata_format) {
        error("ADATA_MERGE only outputs the AnnData object; anndata_format must not be empty")
    }
    template('merge.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def stub_cmd = anndata_format == 'zarr'
        ? "mkdir ${prefix}_outer.zarr ${prefix}_integrate.zarr && touch ${prefix}_outer.zarr/zarr.json ${prefix}_integrate.zarr/zarr.json"
        : "touch ${prefix}_outer.h5ad ${prefix}_integrate.h5ad"
    """
    ${stub_cmd}
    touch versions.yml
    """
}
