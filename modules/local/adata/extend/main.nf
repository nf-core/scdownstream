process ADATA_EXTEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/74/74c1dc16cd8944bf1260dd8420da6ae08cfe0f15ac4d382d24e4857d6b5b616a/data'
        : 'community.wave.seqera.io/library/anndata_pandas_python_pyyaml:f3cd364e43f4f891'}"

    input:
    tuple val(meta), path(base), path(obs, stageAs: 'obs/'), path(var, stageAs: 'var/'), path(obsm, stageAs: 'obsm/'), path(obsp, stageAs: 'obsp/'), path(uns, stageAs: 'uns/'), path(layers, stageAs: 'layers/')
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    tuple val(meta), path("*.csv"), emit: metadata
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    template('extend.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch ${prefix}_metadata.csv
    touch versions.yml
    """
}
