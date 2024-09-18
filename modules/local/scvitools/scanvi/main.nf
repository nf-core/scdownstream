process SCVITOOLS_SCANVI {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'docker.io/nicotru/scvitools-gpu:cuda-12' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata_scvi-tools:54d2eb2f946e0e43':
        'community.wave.seqera.io/library/anndata_scvi-tools:fa9451a13918eae0' }"

    input:
    tuple val(meta), path(h5ad)
    tuple val(meta2), path(reference_model)

    output:
    tuple val(meta), path("*.h5ad") , emit: h5ad
    tuple val(meta), path("*_model"), emit: model
    path "${prefix}.pkl"            , emit: obs
    path "X_${prefix}.pkl"          , emit: obsm
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'scanvi.py'
}
