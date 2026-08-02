process SCVITOOLS_SCVI {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.use_gpu
        ? 'ghcr.io/scverse/scvi-tools:py3.13-cu12-1.4.3-'
        : workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/df/dfb4b54fd5cb5c5624d947914f5ab8ac86eeaa5f762ebc52c034fbd36cf30250/data'
            : 'community.wave.seqera.io/library/scvi-tools:1.4.3--cce8c95b58ececa6'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
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

    output:
    tuple val(meta), path("${prefix}.h5ad")          , emit: h5ad
    tuple val(meta), path("${prefix}_model/model.pt"), emit: model
    path "X_${prefix}.pkl"                           , emit: obsm
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

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
