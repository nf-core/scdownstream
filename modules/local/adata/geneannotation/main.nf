process ADATA_GENEANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/ee30d9b6847a686cb030b33ee3d2a1a82446bed2893bad6c47be1eb8e020f7a1/data'
        : 'community.wave.seqera.io/library/mygene_anndata_python_pyarrow_pyyaml:b60d3ce6f3e9c6aa'}"

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
