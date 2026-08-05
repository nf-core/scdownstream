process SCVITOOLS_SCANVI {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/21/21c048694d2a9bc03c17a310c79c5691d2fc83df62e5e2654791a2105da8b09f/data'
:         'community.wave.seqera.io/library/scvi-tools_pyarrow:1734dedd3c3d134b' }"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    tuple val(meta2), path(reference_model, stageAs: 'reference_model/model.pt')
    val reference_model_type
    tuple val(label_col), val(unlabeled_category)
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
    path "${prefix}.parquet"                             , emit: obs
    path "X_${prefix}.parquet"                           , emit: obsm
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    if ("${h5ad}" == "${prefix}.h5ad") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    template('scanvi.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    mkdir -p ${prefix}_model
    touch ${prefix}_model/model.pt
    touch ${prefix}.parquet
    touch X_${prefix}.parquet
    touch versions.yml
    """
}
