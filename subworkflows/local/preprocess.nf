include { ADATA_READRDS                } from '../../modules/local/adata/readrds'
include { ADATA_READCSV                } from '../../modules/local/adata/readcsv'
include { ADATA_UNIFY                  } from '../../modules/local/adata/unify'
include { SCANPY_PLOTQC as QC_RAW      } from '../../modules/local/scanpy/plotqc'
include { CELDA_DECONTX                } from '../../modules/local/celda/decontx'
include { DOUBLET_DETECTION            } from './doublet_detection'
include { SCANPY_PLOTQC as QC_FILTERED } from '../../modules/local/scanpy/plotqc'

workflow PREPROCESS {

    take:
    ch_samples // channel: [ val(meta), file ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samples = ch_samples.map { meta, file -> [meta, file, file.extension.toLowerCase()] }
        .branch { meta, file, ext ->
            h5ad: ext == "h5ad"
                return [meta, file]
            rds: ext == "rds"
                return [meta, file]
            csv: ext == "csv"
                return [meta, file]
        }

    ADATA_READRDS(ch_samples.rds)
    ch_h5ad = ch_samples.h5ad.mix(ADATA_READRDS.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READRDS.out.versions)

    ADATA_READCSV(ch_samples.csv)
    ch_h5ad = ch_samples.h5ad.mix(ADATA_READCSV.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READCSV.out.versions)

    ADATA_UNIFY(ch_h5ad)
    ch_h5ad = ADATA_UNIFY.out.h5ad
    ch_versions = ch_versions.mix(ADATA_UNIFY.out.versions)

    QC_RAW(ch_h5ad)
    ch_multiqc_files = ch_multiqc_files.mix(QC_RAW.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_RAW.out.versions)

    // Ambient RNA removal
    if (params.ambient_removal == 'decontx') {
        CELDA_DECONTX(ch_h5ad)
        ch_h5ad = CELDA_DECONTX.out.h5ad
        ch_versions = CELDA_DECONTX.out.versions
    }

    DOUBLET_DETECTION(ch_h5ad)
    ch_h5ad = DOUBLET_DETECTION.out.h5ad
    ch_multiqc_files = ch_multiqc_files.mix(DOUBLET_DETECTION.out.multiqc_files)
    ch_versions = ch_versions.mix(DOUBLET_DETECTION.out.versions)

    QC_FILTERED(ch_h5ad)
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    emit:
    h5ad          = ch_h5ad

    multiqc_files = ch_multiqc_files
    versions      = ch_versions                     // channel: [ versions.yml ]
}
