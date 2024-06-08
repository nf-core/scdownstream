process ADATA_READRDS {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:65d6b87e3a6a9270':
        'community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:693764013a2a63f5' }"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'readrds.py'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        adata: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
