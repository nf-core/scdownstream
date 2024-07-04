process CELLTYPES_SINGLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-celldex_bioconductor-singlecellexperiment_bioconductor-singler_anndata:d0dfcaede2417581':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-celldex_bioconductor-singlecellexperiment_bioconductor-singler_anndata:d6a21ee363999d21' }"

    input:
    tuple val(meta), path(h5ad)
    val(reference)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path("*.pkl")                  , emit: obs
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'singleR.py'
}
