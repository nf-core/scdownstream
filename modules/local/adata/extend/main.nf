process ADATA_EXTEND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/68/68042761ce4b44cb86292bc2f733dfc0699a67fc3d4531761ddc1c921ad3b485/data'
:         'community.wave.seqera.io/library/anndata_pandas_pyarrow_python_pyyaml:684cdfb0284a3431' }"

    input:
    tuple val(meta), path(base), path(obs, stageAs: 'obs/'), path(var, stageAs: 'var/'), path(obsm, stageAs: 'obsm/'), path(obsp, stageAs: 'obsp/'), path(uns, stageAs: 'uns/'), path(layers, stageAs: 'layers/')

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    tuple val(meta), path("*.csv"), emit: metadata
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('extend.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch ${prefix}_metadata.csv
    touch versions.yml
    """
}
