process SCVITOOLS_SCVI {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${task.accelerator
        ? 'ghcr.io/scverse/scvi-tools:py3.13-cu12-1.4.3-'
        : workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9ea53a1a670cb5d83abab8b885034e786fed9ff52cb683aa3f1cc9bd5aae2d93/data'
            : 'community.wave.seqera.io/library/scvi-tools:1.4.3--05da7945b26925f8'}"

    input:
    tuple val(meta), path(anndata, arity: 1)
    tuple val(meta2), path(reference_model, stageAs: 'reference_model/model.pt')
    val batch_col
    val categorical_covariates
    val continuous_covariates
    val n_hidden
    val n_layers
    val n_latent
    val dispersion
    val gene_likelihood
    val max_epochs
    val use_observed_lib_size
    val anndata_format

    output:
    tuple val(meta), path("${prefix}.{h5ad,zarr}"), emit: anndata, optional: true
    tuple val(meta), path("${prefix}_model/model.pt"), emit: model
    path "X_${prefix}.pkl", emit: obsm
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!(anndata_format in ['h5ad', 'zarr', ''])) {
        error("anndata_format must be 'h5ad', 'zarr' or '', got '${anndata_format}'")
    }
    if (anndata_format && "${prefix}.${anndata_format}" == "${anndata}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }
    template('scvi.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def anndata_stub = anndata_format == 'zarr'
        ? "mkdir ${prefix}.zarr && touch ${prefix}.zarr/zarr.json"
        : anndata_format == 'h5ad' ? "touch ${prefix}.h5ad" : ''
    """
    ${anndata_stub}
    mkdir -p ${prefix}_model
    touch ${prefix}_model/model.pt
    touch X_${prefix}.pkl
    touch versions.yml
    """
}
