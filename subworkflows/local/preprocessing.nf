workflow PREPROCESSING {

    take:
    ch_datasets // channel: [ val(meta), file ]

    main:

    ch_versions = Channel.empty()

    ch_datasets = ch_datasets.map { meta, file -> [meta, file, file.extension.toLowerCase()] }
        .branch { meta, file, ext ->
            h5ad: ext == "h5ad"
                return [meta, file]
            rds: ext == "rds"
                return [meta, file]
        }

    ch_datasets.h5ad.view()

    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}

