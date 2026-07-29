process LIANA_BYSAMPLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/00a0832cf7999e945fe1eea426fb5a59d4718af37c8e97b1768d7257904bfadb/data' :
        'community.wave.seqera.io/library/liana_pyyaml:c2bc77d30223be1f' }"

    input:
    tuple val(meta), path(h5ad)
    val obs_key
    val context_key
    val n_perms
    val max_cells
    val subsample_strategy
    val subsample_seed

    output:
    tuple val(meta), path("*.csv.gz")      , emit: results, optional: true
    tuple val(meta), path("*_contexts.tsv"), emit: contexts, optional: true
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'bysample.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip > "${prefix}.csv.gz"
    touch "${prefix}_contexts.tsv"
    touch "versions.yml"
    """
}
