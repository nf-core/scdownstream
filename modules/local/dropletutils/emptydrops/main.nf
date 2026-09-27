process DROPLETUTILS_EMPTYDROPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24bb1be29bfea3307665b3f1f1d294e3a35f871ee55f1728d4ec63c585fa755a/data'
        : 'community.wave.seqera.io/library/bioconductor-anndatar_bioconductor-dropletutils_bioconductor-rhdf5_bioconductor-singlecellexperiment:1a7241b7686e28be'}"

    input:
    tuple val(meta), path(h5ad)
    val lower
    val fdr

    output:
    tuple val(meta), path("${prefix}_barcodes.csv"), emit: barcodes
    tuple val(meta), path("${prefix}_results.csv"), emit: results
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('emptydrops.R')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_barcodes.csv
    touch ${prefix}_results.csv
    touch versions.yml
    """
}
