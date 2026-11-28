process ADATA_GENEANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/da/dadc22a2e2add94e0cbcb29e01199de9f71594c4eaa1f6286e693a64139900d7/data'
        : 'community.wave.seqera.io/library/mygene_anndata_python_pyarrow_pyyaml:00f473d9060bc114'}"

    input:
    tuple val(meta), path(h5ad)
    path gene_table
    val species

    output:
    tuple val(meta), path("*.parquet"), emit: var
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_gene_annotation"
    template('geneannotation.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_gene_annotation"
    """
    echo ",gene_biotype,chromosome,biotype_group,gene_annotation_ambiguous" > ${prefix}.csv
    touch ${prefix}.parquet
    touch versions.yml
    """
}
