process ADATA_EXTEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0d347e0d943cb8d0624a8c73247db0c8e96e7c130fe0fb3125aeb7d1b140e79/data' :
        'community.wave.seqera.io/library/anndata_python_pyyaml:5ac0fd4d280528a5' }"

    input:
    tuple (
        val(meta),
        path(base),
        path(obs, stageAs: 'obs/'),
        path(var, stageAs: 'var/'),
        path(obsm, stageAs: 'obsm/'),
        path(obsp, stageAs: 'obsp/'),
        path(uns, stageAs: 'uns/'),
        path(layers, stageAs: 'layers/')
    )

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.csv") , emit: metadata
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('extend.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_metadata.csv
    touch versions.yml
    """
}
