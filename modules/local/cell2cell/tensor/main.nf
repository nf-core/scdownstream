process CELL2CELL_TENSOR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/10d2843fde28026d17f6faffc136aaa352f85d91fec25d658f793e062f14f349/data' :
        'community.wave.seqera.io/library/liana_cell2cell_pyyaml:50d977636ad3e43e' }"

    input:
    tuple val(meta), path(liana_bysample, stageAs: 'liana_bysample.csv.gz'), path(contexts, stageAs: 'contexts.tsv')
    val context_key
    val rank
    val seed
    val species

    output:
    tuple val(meta), path("*.png")  , emit: plots, optional: true
    path("*_loadings_*.csv")        , emit: loadings, optional: true
    path("*_pathway_enrichment.csv"), emit: enrichment, optional: true
    path("*_tensor.pkl")            , emit: tensor, optional: true
    path "*_mqc.json"               , emit: multiqc_files, optional: true
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix      = task.ext.prefix ?: "${meta.id}"
    integration = meta.integration ?: 'integration'
    template 'tensor_c2c.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_factor_1_loadings_product.png"
    touch "${prefix}_tensor_factors.png"
    touch "${prefix}_loadings_lr_clustermap.png"
    touch "${prefix}_loadings_contexts_clustermap.png"
    touch "${prefix}_pathway_enrichment_dotplot.png"
    touch "${prefix}_loadings_contexts.csv"
    touch "${prefix}_loadings_ligand_receptor_pairs.csv"
    touch "${prefix}_loadings_sender_cells.csv"
    touch "${prefix}_loadings_receiver_cells.csv"
    touch "${prefix}_tensor.pkl"
    echo '{}' > "${prefix}_factor_1_mqc.json"
    echo '{}' > "${prefix}_tensor_factors_mqc.json"
    echo '{}' > "${prefix}_loadings_lr_clustermap_mqc.json"
    echo '{}' > "${prefix}_loadings_contexts_clustermap_mqc.json"
    echo '{}' > "${prefix}_pathway_enrichment_dotplot_mqc.json"
    touch "versions.yml"
    """
}
