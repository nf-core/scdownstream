process CELLBENDER_REMOVEBACKGROUND {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_long'
    label 'process_gpu'

    conda "${moduleDir}/${task.accelerator ? 'environment.gpu.yml' : 'environment.yml'}"
    container "${task.accelerator
        ? 'community.wave.seqera.io/library/cellbender_python_pytorch_cuda-version_webcolors:c8b51e8dfdbcd308'
        : workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5b5c128de14d5f9ae914eeab43aea7aee2a0994f3a57fd98e1a1b6656df34d34/data'
            : 'community.wave.seqera.io/library/cellbender_python_pytorch_webcolors:4643692df2dc7df0'}"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("${prefix}.h5"), emit: h5
    tuple val(meta), path("${prefix}_filtered.h5"), emit: filtered_h5
    tuple val(meta), path("${prefix}_posterior.h5"), emit: posterior_h5
    tuple val(meta), path("${prefix}_cell_barcodes.csv"), emit: barcodes
    tuple val(meta), path("${prefix}_metrics.csv"), emit: metrics
    tuple val(meta), path("${prefix}_report.html"), emit: report, optional: true
    tuple val(meta), path("${prefix}.pdf"), emit: pdf
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val(meta), path("ckpt.tar.gz"), emit: checkpoint
    // The Bioconda 0.4.0 build retains the 0.3.2 Python distribution metadata.
    tuple val("${task.process}"), val('cellbender'), eval('python3 -c \'from importlib.metadata import version; from pathlib import Path; import sys; print(next((record.stem.split("-")[1] for record in Path(sys.prefix, "conda-meta").glob("cellbender-*.json")), version("cellbender")))\''), emit: versions_cellbender, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args = task.ext.args ?: ""
    use_gpu = task.accelerator ? "--cuda" : ""
    """
        TMPDIR=. cellbender remove-background \
            ${args} \
            --cpu-threads ${task.cpus} \
            ${use_gpu} \
            --input ${h5ad} \
            --output ${prefix}.h5
        """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
        touch "${prefix}.h5"
        touch "${prefix}_filtered.h5"
        touch "${prefix}_posterior.h5"
        touch "${prefix}_cell_barcodes.csv"
        touch "${prefix}_metrics.csv"
        touch "${prefix}_report.html"
        touch "${prefix}.pdf"
        touch "${prefix}.log"
        echo "" | gzip > ckpt.tar.gz
        """
}
