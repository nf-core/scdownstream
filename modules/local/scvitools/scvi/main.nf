process SCVITOOLS_SCVI {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.use_gpu
        ? 'ghcr.io/scverse/scvi-tools:py3.12-cu12-1.3.3-'
        : workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8764e4208e9639a54d636fc65c839c55dedbfd68def57baea90d1d2007d6a7f/data'
            : 'community.wave.seqera.io/library/scvi-tools:1.3.3--df115aabdccb7d6b'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    tuple val(meta2), path(reference_model, stageAs: 'reference_model/model.pt')
    val batch_col
    val categorical_covariates
    val continuous_covariates

    output:
    tuple val(meta), path("${prefix}.h5ad")          , emit: h5ad
    tuple val(meta), path("${prefix}_model/model.pt"), emit: model
    path "X_${prefix}.pkl"                           , emit: obsm
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    n_hidden = task.ext.n_hidden ?: 128
    n_layers = task.ext.n_layers ?: 2
    n_latent = task.ext.n_latent ?: 30
    dispersion = task.ext.dispersion ?: 'gene'
    gene_likelihood = task.ext.gene_likelihood ?: 'zinb'
    max_epochs = task.ext.max_epochs ?: null

    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    template('scvi.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    mkdir -p ${prefix}_model
    touch ${prefix}_model/model.pt
    touch X_${prefix}.pkl
    touch versions.yml
    """
}
