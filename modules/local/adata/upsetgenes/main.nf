process ADATA_UPSETGENES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scanpy_upsetplot:7f44fbc2599e351d':
        'community.wave.seqera.io/library/scanpy_upsetplot:6108c2e3903b49e2' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.png"), emit: plot
    path("*_mqc.json")            , emit: multiqc_files
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix    ?: "${meta.id}"
    split_col = task.ext.split_col ?: 'sample'
    template 'upsetplot.py'
}
