process SCANPY_RANKGENESGROUPS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a497dfbc392812dd3fec43e8c3cdc653e62f868865c25be8805f38ec908ecc69/data'
:         'community.wave.seqera.io/library/adjusttext_anndata_pandas_python_pruned:a229cbba149bfd01' }"

    input:
    tuple val(meta), path(h5ad)
    val(obs_key)
    tuple val(filter_col), val(filter_val)
    val(method)
    val(rank_key)
    path interesting_genes

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad, optional: true
    path "*.pkl"                   , emit: uns, optional: true
    path "*.png"                   , emit: plots, optional: true
    path "*.csv"                   , emit: markers, optional: true
    path "*_mqc.json"              , emit: multiqc_files, optional: true
    path "versions.yml"            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('rank_genes_groups.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "${prefix}.pkl"
    touch "${prefix}.png"
    touch "${prefix}_dotplot.png"
    touch "${prefix}_volcano.png"
    touch "${prefix}_mqc.json"
    touch "${prefix}_dotplot_mqc.json"
    touch "${prefix}_volcano_mqc.json"
    touch "versions.yml"
    """
}
