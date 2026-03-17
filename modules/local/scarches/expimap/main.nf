process SCARCHES_EXPIMAP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
            ? 'https://wave.seqera.io/view/builds/bd-c86ba663882352b6_1'
            : 'community.wave.seqera.io/library/pandas_python_pip_anndata_pruned:c86ba663882352b6'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)
    tuple val(meta2), path(reference_model)
    val(condition_col)
    val(counts_layer)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "X_${prefix}.pkl", emit: obsm
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    hidden_layer_sizes = task.ext.hidden_layer_sizes ?: [256, 256, 256]
    recon_loss = task.ext.recon_loss ?: "nb"
    n_epochs = task.ext.n_epochs ?: 400
    alpha_epoch_anneal = task.ext.alpha_epoch_anneal ?: 100
    alpha = task.ext.alpha ?: 0.7
    alpha_kl = task.ext.alpha_kl ?: 0.5
    use_early_stopping = task.ext.use_early_stopping ?: true
    if ("${prefix}.h5ad" == "${h5ad}") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    template('expimap.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch X_${prefix}.pkl
    touch versions.yml
    """
}
