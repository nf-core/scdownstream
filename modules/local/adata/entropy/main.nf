process ADATA_ENTROPY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/38/386e532a346d91ce6e59118da064b9367024f807fb7a52d1e13457a0d01107f2/data'
        : 'community.wave.seqera.io/library/python_pyyaml_scanpy:ab265932a4ff2ebf'}"

    input:
    tuple val(meta), path(anndata)
    val group_col
    val entropy_col
    val plot_basis
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    path "${prefix}.pkl", emit: obs
    path "${prefix}.png", emit: plots, optional: true
    path "${prefix}_mqc.json", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    template('entropy.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_entropy"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    touch ${prefix}.pkl

    if [ ${plot_basis ? 'true' : 'false'} ]; then
        touch ${prefix}.png
        touch ${prefix}_mqc.json
    fi

    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
