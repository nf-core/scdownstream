process CUSTOM_VOLCANOPLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/2367404f0aa4c23d7e0c15e24ccc53444bb650fbae12282a680db6c6ad3bdfbb/data':
        'community.wave.seqera.io/library/custom_volcanoplot:2e9b936eb387e94b' }"

    input:
    tuple val(meta), path(parquet)
    path interesting_genes

    output:
    path "*.png"       , emit: plots        , optional: true
    path "*_mqc.json"  , emit: multiqc_files, optional: true
    path "versions.yml", emit: versions     , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'volcanoplot.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def plot_stem = parquet ? parquet.baseName.replaceFirst('_results$', '') : prefix
    """
    touch ${plot_stem}_volcano.png
    touch ${plot_stem}_volcano_mqc.json
    touch versions.yml
    """
}
