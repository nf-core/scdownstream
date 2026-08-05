process LIANA_RANKAGGREGATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/35fa059b36b97f47fb2dafbb573ab9505647f08dd673efbd9913f74f1556ce68/data'
:         'community.wave.seqera.io/library/liana_pyyaml_pyarrow:5fc198487f068dd6' }"

    input:
    tuple val(meta), path(h5ad)
    val n_perms
    val max_cells
    val subsample_strategy
    val subsample_seed

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad         , optional: true
    path("*.parquet")                  , emit: uns          , optional: true
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
    touch "${prefix}.parquet"
    touch "${prefix}_dotplot.png"
    touch "${prefix}_circle.png"
    touch "${prefix}_tileplot.png"
    touch "${prefix}_dotplot_mqc.json"
    touch "${prefix}_circle_mqc.json"
    touch "${prefix}_tileplot_mqc.json"
    touch "versions.yml"
    """
}
