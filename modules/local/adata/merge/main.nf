process ADATA_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8e2e626f2e43e75b0f218a91d16c13ac7efd15c6e0618dbc8ce06cd84218ca99/data'
        : 'community.wave.seqera.io/library/python-igraph_python_pyyaml_scanpy:1ea7adbfe9ea8cad'}"

    input:
    tuple val(meta), path(h5ads, stageAs: 'input/sample_?.h5ad')
    tuple val(meta2), path(base)
    val force_obs_cols

    output:
    tuple val(meta), path("*_outer.h5ad"), emit: outer
    tuple val(meta), path("*_integrate.h5ad"), emit: integrate
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('merge.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_outer.h5ad
    touch ${prefix}_integrate.h5ad
    touch versions.yml
    """
}
