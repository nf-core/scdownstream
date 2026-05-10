include { SCANPY_READH5 } from '../../../modules/local/scanpy/readh5'
include { ADATA_READRDS } from '../../../modules/local/adata/readrds'
include { ADATA_READCSV } from '../../../modules/local/adata/readcsv'

workflow LOAD_H5AD {
    take:
    ch_samples // channel: [ meta, {h5ad/rds}, {h5ad/rds} ]

    main:
    ch_versions = channel.empty()
    ch_files = channel.empty()
    ch_h5ad = channel.empty()

    ch_files = ch_files.mix(ch_samples
        .map {
            meta, filtered, _unfiltered
            -> [meta + [type: 'filtered'], filtered] }
        .filter { _meta, filtered -> filtered }
    )
    ch_files = ch_files.mix(ch_samples
        .map {
            meta, _filtered, unfiltered
            -> [meta + [type: 'unfiltered'], unfiltered] }
        .filter { _meta, unfiltered -> unfiltered }
    )
    ch_files = ch_files
        .map { meta, file -> [meta, file, file.extension.toLowerCase()] }
        .branch { meta, file, ext ->
            h5ad: ext == "h5ad"
            return [meta, file]
            h5: ext == "h5"
            return [meta, file]
            rds: ext == "rds"
            return [meta, file]
            csv: ext == "csv"
            return [meta, file]
        }

    ch_h5ad = ch_h5ad.mix(ch_files.h5ad)

    SCANPY_READH5 (
        ch_files.h5
    )
    ch_h5ad = ch_h5ad.mix(SCANPY_READH5.out.h5ad)
    ch_versions = ch_versions.mix(SCANPY_READH5.out.versions)

    ADATA_READRDS (
        ch_files.rds
    )
    ch_h5ad = ch_h5ad.mix(ADATA_READRDS.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READRDS.out.versions)

    ADATA_READCSV (
        ch_files.csv
    )
    ch_h5ad = ch_h5ad.mix(ADATA_READCSV.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READCSV.out.versions)

    ch_output = ch_samples
        .map { meta, _filtered, _unfiltered -> [meta.id, meta] }
        .join(ch_h5ad
                .filter { meta, _h5ad -> meta.type == 'filtered' }
                .map { meta, filtered -> [meta.id, filtered] },
            failOnMismatch: false,
            remainder: true,
        )
        .join(ch_h5ad
                .filter { meta, _h5ad -> meta.type == 'unfiltered' }
                .map { meta, unfiltered -> [meta.id, unfiltered] },
            failOnMismatch: false,
            remainder: true,
        )
        .map {
            _id, meta, filtered, unfiltered
            -> [meta, filtered ?: [], unfiltered ?: []]
        }

    emit:
    h5ad     = ch_output   // channel: [ meta, h5ad, h5ad ]
    versions = ch_versions // channel: [ versions.yml ]
}
