process SCVITOOLS_SCVI {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/${task.accelerator ? 'environment.gpu.yml' : 'environment.yml'}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/0072fc1c3e82f329a0f92c8544e6f9b2129ff08f7cbfa32b50f8c3df1b54f3f3/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba88d310e8f5d69d321b94bb73012e8f5571ec4aee16e1047e7dd378b606fb2a/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/python_scvi-tools_pyarrow_pytorch-gpu_cuda-version@sha256:7bd8b3dd714d9e419a67bbbaab1d92fa5823fd370b307c7e67d38a02ed36cbb1'
            : 'community.wave.seqera.io/library/pyarrow_scvi-tools:002a45f506b3d900')}"

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
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    tuple val(meta), path("${prefix}_model/model.pt"), emit: model
    path "X_${prefix}.parquet", emit: obsm
    path "versions.yml", emit: versions, topic: versions

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
    touch X_${prefix}.parquet
    touch versions.yml
    """
}
