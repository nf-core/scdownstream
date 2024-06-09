process CELDA_DECONTX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-celda_anndata:de20d5cbd4f86aa6':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-celda_anndata:31bbf686a87fe0aa' }"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    batch_col = task.ext.batch_col ?: "batch"
    template 'decontx.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        celda: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
