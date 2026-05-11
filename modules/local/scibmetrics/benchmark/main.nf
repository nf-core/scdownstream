process SCIBMETRICS_BENCHMARK {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/342afc015653aee4762076abcb281f82d815e948644bbb05a1af688eb1980506/data'
        : 'community.wave.seqera.io/library/python_pyyaml_pip_scib-metrics:26d5406663c7db7d'}"

    input:
    tuple val(meta), path(h5ad, arity: 1)

    output:
    path "${prefix}_metrics.tsv", emit: metrics
    path "${prefix}_mqc.json"   , emit: multiqc_files
    path "versions.yml"         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ''
    template('benchmark.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metrics.tsv
    touch ${prefix}_mqc.json
    touch versions.yml
    """
}
