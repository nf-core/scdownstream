process ADATA_MYGENE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/ccc5af513e241deec87c61698bd3d49f458a78938474ff6e2b9a006dfc51690a/data'
        : 'community.wave.seqera.io/library/mygene_anndata_python_pyyaml:0b27d9b10dbf88bf'}"

    input:
    tuple val(meta), path(anndata)
    val species
    val input_col
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
        error("ADATA_MYGENE only outputs the AnnData object; anndata_format must not be empty")
    }
    if (anndata_format && "${prefix}.${anndata_format}" == "${anndata}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('mygene.py')

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
