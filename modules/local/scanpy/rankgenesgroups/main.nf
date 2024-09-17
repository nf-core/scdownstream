process SCANPY_RANKGENESGROUPS {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ task.ext.use_gpu ? 'ghcr.io/scverse/rapids_singlecell:v0.10.8' :
        workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scanpy:1.10.1--ea08051addf267ac':
        'community.wave.seqera.io/library/scanpy:1.10.1--0c8c97148fc05558' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "*.pkl"                   , emit: uns
    path "*.png"                   , emit: plots
    path "*_mqc.json"              , emit: multiqc_files
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    obs_key = task.ext.obs_key ?: "leiden"
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'rank_genes_groups.py'
}
