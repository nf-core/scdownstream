process DOUBLET_REMOVAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata_upsetplot:ca171abd2131e1ba':
        'community.wave.seqera.io/library/anndata_upsetplot:59615c04a744efcb' }"

    input:
    tuple val(meta), path(h5ad)
    path(predictions)
    val(threshold)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path("*_mqc.png")              , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix    ?: "${meta.id}"
    template 'doublet_removal.py'
}
