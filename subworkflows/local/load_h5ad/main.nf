include { SCANPY_READH5 } from '../../../modules/local/scanpy/readh5'
include { ADATA_READRDS } from '../../../modules/local/adata/readrds'
include { ADATA_READCSV } from '../../../modules/local/adata/readcsv'

def anndataFormat(meta, h5ad_required) {
    (meta.type == 'unfiltered' || h5ad_required) ? 'h5ad' : 'zarr'
}

workflow LOAD_H5AD {
    take:
    ch_samples // channel: [ meta, {h5ad/rds}, {h5ad/rds} ]
    h5ad_required //   value: boolean, true if an h5ad-only consumer (ambient correction, hugo-unifier) reads the filtered objects

    main:
    ch_files = channel.empty()
    ch_h5ad = channel.empty()

    ch_files = ch_files.mix(
        ch_samples.map { meta, filtered, _unfiltered ->
            [meta + [type: 'filtered'], filtered]
        }.filter { _meta, filtered -> filtered }
    )
    ch_files = ch_files.mix(
        ch_samples.map { meta, _filtered, unfiltered ->
            [meta + [type: 'unfiltered'], unfiltered]
        }.filter { _meta, unfiltered -> unfiltered }
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

    ch_h5 = ch_files.h5.multiMap { meta, file ->
        file: [meta, file]
        anndata_format: anndataFormat(meta, h5ad_required)
    }
    SCANPY_READH5(
        ch_h5.file,
        ch_h5.anndata_format,
    )
    ch_h5ad = ch_h5ad.mix(SCANPY_READH5.out.anndata)

    ch_rds = ch_files.rds.multiMap { meta, file ->
        file: [meta, file]
        anndata_format: anndataFormat(meta, h5ad_required)
    }
    ADATA_READRDS(
        ch_rds.file,
        ch_rds.anndata_format,
    )
    ch_h5ad = ch_h5ad.mix(ADATA_READRDS.out.anndata)

    ch_csv = ch_files.csv.multiMap { meta, file ->
        file: [meta, file]
        anndata_format: anndataFormat(meta, h5ad_required)
    }
    ADATA_READCSV(
        ch_csv.file,
        ch_csv.anndata_format,
    )
    ch_h5ad = ch_h5ad.mix(ADATA_READCSV.out.anndata)

    ch_output = ch_samples
        .map { meta, _filtered, _unfiltered -> [meta.id, meta] }
        .join(
            ch_h5ad.filter { meta, _h5ad -> meta.type == 'filtered' }.map { meta, filtered -> [meta.id, filtered] },
            failOnMismatch: false,
            remainder: true,
        )
        .join(
            ch_h5ad.filter { meta, _h5ad -> meta.type == 'unfiltered' }.map { meta, unfiltered -> [meta.id, unfiltered] },
            failOnMismatch: false,
            remainder: true,
        )
        .map { _id, meta, filtered, unfiltered ->
            [meta, filtered ?: [], unfiltered ?: []]
        }

    emit:
    anndata = ch_output // channel: [ meta, anndata, anndata ]
}
