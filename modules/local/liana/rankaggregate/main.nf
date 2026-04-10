process LIANA_RANKAGGREGATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e83ce3d883af5a7be7f05e740fffaf0fff63b5f13e9bf175af9465e91c8cfda2/data' :
        'community.wave.seqera.io/library/liana_pyyaml:776fdd7103df146d' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path("*.pkl")                  , emit: uns, optional: true
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args    = task.ext.args   ?: ''
    obs_key = meta.obs_key ?: "leiden"
    prefix  = task.ext.prefix ?: "${meta.id}"
    template 'rank_aggregate.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "versions.yml"
    """
}
