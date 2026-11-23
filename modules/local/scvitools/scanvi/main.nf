process SCVITOOLS_SCANVI {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/${task.accelerator ? 'environment.gpu.yml' : 'environment.yml'}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/07/0703860e978cf8a1f420cd6f39a5d32a3e3bb255fbf547ae47f673ebc5b76fbe/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/4867a52f4a6f4a5216a3c21e2aead3f7eee65495b517c7b332e3bdf8cea13b30/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/python_scvi-tools_pyarrow_pytorch-gpu_cuda-version@sha256:069836e30ad5879e2ca3d2fc1c0904fdbff649add0cd19ac3c31afddc55aee4e'
            : 'community.wave.seqera.io/library/python_scvi-tools_pyarrow_pytorch-cpu:f0c65d8a2e3eacb6')}"

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
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}_model/model.pt"), emit: model
    path "${prefix}.parquet", emit: obs
    path "X_${prefix}.parquet", emit: obsm
    path "versions.yml", emit: versions, topic: versions

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
