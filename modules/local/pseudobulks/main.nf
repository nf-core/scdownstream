process ADATA_PSEUDOBULKS {
    tag "$meta.id"
    label 'process_medium'
    //label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/scanpy:1.10.4--c2d474f46255931c':
        'community.wave.seqera.io/library/scanpy:1.10.4--f905699eb17b6536' }"

    input:
    tuple val(meta),  path(h5ad)
    val(groups)
    val(mode)
    val(replicates_per_sample)
    val(min_cells)


    output:
    tuple val(meta), path("*.h5ad")    , emit: h5ad
    path "pseudobulks.pkl"             , emit: pseudobulks
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'pseudobulk.py'
}
