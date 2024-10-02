process SCANPY_PAGA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python-igraph_scanpy:a7114d55b0af3893':
        'community.wave.seqera.io/library/python-igraph_scanpy:5f677450e42211ef' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path("*.pkl")                  , emit: uns
    path("*.npy")                  , emit: obsp
    path("*.png")                  , emit: plot
    path("*_mqc.json")             , emit: multiqc_files
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    obs_key = meta.obs_key ?: "leiden"
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'paga.py'
}
