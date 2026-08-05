process ADATA_MERGEEMBEDDINGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f4cb85cb4a864e74997f047afb5f16ecd6e106cc8793edfc3e91945501b41bbd/data'
:         'community.wave.seqera.io/library/python_pyyaml_scanpy_pyarrow:4be15c4d717947be' }"

    input:
    tuple val(meta), val(integration_key), path(integrated, stageAs: 'integrated.h5ad'), path(base), path(combined)

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path ("${prefix}.parquet")                 , emit: obs, optional: true
    path ("X_${prefix}.parquet")               , emit: obsm
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('merge_embeddings.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}.parquet
    touch X_${prefix}.parquet
    touch versions.yml
    """
}
