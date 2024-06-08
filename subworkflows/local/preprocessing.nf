include { ADATA_READRDS                } from '../../modules/local/adata/readrds'
include { SCANPY_PLOTQC as QC_RAW      } from '../../modules/local/scanpy/plotqc'
include { SCANPY_PLOTQC as QC_FILTERED } from '../../modules/local/scanpy/plotqc'

workflow PREPROCESSING {

    take:
    ch_datasets // channel: [ val(meta), file ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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

    QC_RAW(ch_h5ad)
    ch_multiqc_files = ch_multiqc_files.mix(QC_RAW.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_RAW.out.versions)

    QC_FILTERED(ch_h5ad)
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    emit:

    multiqc_files = ch_multiqc_files
    versions      = ch_versions                     // channel: [ versions.yml ]
}
