process LIANA_RANKAGGREGATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00a0832cf7999e945fe1eea426fb5a59d4718af37c8e97b1768d7257904bfadb/data'
:         'community.wave.seqera.io/library/liana_pyyaml:c2bc77d30223be1f' }"

    input:
    tuple val(meta), path(h5ad)
    val n_perms
    val max_cells
    val subsample_strategy
    val subsample_seed

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad         , optional: true
    path("*.pkl")                  , emit: uns          , optional: true
    path("*.png")                  , emit: plots        , optional: true
    path("*_mqc.json")             , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    obs_key = meta.obs_key ?: "leiden"
    prefix  = task.ext.prefix ?: "${meta.id}"
    template 'rank_aggregate.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "${prefix}_dotplot.png"
    touch "${prefix}_circle.png"
    touch "${prefix}_tileplot.png"
    touch "${prefix}_dotplot_mqc.json"
    touch "${prefix}_circle_mqc.json"
    touch "${prefix}_tileplot_mqc.json"
    touch "versions.yml"
    """
}
