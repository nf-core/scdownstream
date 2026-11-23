process CELLBENDER_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b5c128de14d5f9ae914eeab43aea7aee2a0994f3a57fd98e1a1b6656df34d34/data'
        : 'community.wave.seqera.io/library/cellbender_python_pytorch_webcolors:4643692df2dc7df0'}"

    input:
    tuple val(meta), path(filtered), path(unfiltered), path(cellbender_h5)
    val output_layer_name

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: h5ad
    path "versions.yml", emit: versions_cellbender, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    output_layer = output_layer_name ?: "cellbender"

    """
    echo ${output_layer}
    """

    template('merge.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import platform; print(platform.python_version())')
        cellbender: \$(python3 -c 'from importlib.metadata import version; from pathlib import Path; import sys; print(next((record.stem.split("-")[1] for record in Path(sys.prefix, "conda-meta").glob("cellbender-*.json")), version("cellbender")))')
    END_VERSIONS
    """
}
