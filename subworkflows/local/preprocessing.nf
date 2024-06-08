include { ADATA_READRDS } from '../../modules/local/adata/readrds'

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

    ADATA_READRDS(ch_datasets.rds)

    ch_versions = ch_versions.mix(ADATA_READRDS.out.versions)

    ch_h5ad = ch_datasets.h5ad.mix(ADATA_READRDS.out.h5ad)

    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}

