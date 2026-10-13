process SYMPHONY_HARMONYINTEGRATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f801c41ede6dcf6eaadb994993323e9a775676245acadcf25970e963b2549f05/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy_pip_symphonypy:30c480bf2aefe534'}"

    input:
    tuple val(meta), path(anndata)
    val batch_col
    val counts_layer
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    tuple val(meta), path("${prefix}_reference.h5ad"), emit: reference
    path "X_${prefix}.pkl", emit: obsm
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (anndata_format && "${prefix}.${anndata_format}" == "${anndata}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('harmonyintegrate.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch ${prefix}_reference.h5ad
    touch X_${prefix}.pkl
    touch versions.yml
    """
}
