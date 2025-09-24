process SCANPY_RANKGENESGROUPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/baee2c1ee0f6cd0b6a18a6c71bad03370139a77e53cad06464b065f795d52cd0/data'
        : 'community.wave.seqera.io/library/pyyaml_scanpy:a3a797e09552fddc'}"

    input:
    tuple val(meta), path(h5ad), path(cluster_csv)

    output:
    tuple val(meta), path(prefix), emit: outdir
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path "*.pkl", emit: uns, optional: true
    path "*.png", emit: plots, optional: true
    path "*_mqc.json", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: cluster_csv.baseName
    sample_group_col = task.ext.sample_group_col ?: null
    method = task.ext.method ?: 'wilcoxon'
    template('rank_genes_groups.py')
}
