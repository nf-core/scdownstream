process SOUPX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_leidenalg_python-igraph_pruned:c2538291aadd50cb':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_leidenalg_python-igraph_pruned:1b1b2ad4205f41be' }"

    input:
    tuple val(meta), path(h5ad), path(raw)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'soupx.py'
}
