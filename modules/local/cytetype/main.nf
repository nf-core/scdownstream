process CYTETYPE {
    tag "${meta.id}"
    label 'process_medium'

    secret secrets.CYTETYPE_API_KEY ? ['CYTETYPE_API_KEY'] : ''

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4db5f29d0da08e2bdce05993ae80dc04b7d5e0a032d426ca0cfffb495a8a8841/data'
        : 'community.wave.seqera.io/library/anndata_pandas_python_pyyaml_pruned:74d512b309c5b8b5'}"

    input:
    tuple val(meta), path(anndata, stageAs: 'input/*')
    val symbol_col
    val study_context
    val group_key
    val rank_key
    val integration
    val resolution
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    path "${prefix}.pkl", emit: obs
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (anndata_format && "${prefix}.${anndata_format}" == anndata.name) {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('cytetype.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch ${prefix}.pkl
    touch versions.yml
    """
}
