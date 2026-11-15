process ADATA_GENEANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cc/ccc5af513e241deec87c61698bd3d49f458a78938474ff6e2b9a006dfc51690a/data'
:         'community.wave.seqera.io/library/mygene_anndata_python_pyyaml:0b27d9b10dbf88bf' }"

    input:
    tuple val(meta), path(h5ad)
    path(gene_table)
    val(species)

    output:
    tuple val(meta), path("*.pkl"), emit: var
    tuple val(meta), path("*.csv"), emit: csv
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_gene_annotation"
    template 'geneannotation.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_gene_annotation"
    """
    echo ",gene_biotype,chromosome,biotype_group,gene_annotation_ambiguous" > ${prefix}.csv
    touch ${prefix}.pkl
    touch versions.yml
    """
}
