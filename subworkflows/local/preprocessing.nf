include { ADATA_READRDS                } from '../../modules/local/adata/readrds'
include { ADATA_UNIFY                  } from '../../modules/local/adata/unify'
include { SCANPY_PLOTQC as QC_RAW      } from '../../modules/local/scanpy/plotqc'
include { CELDA_DECONTX                } from '../../modules/local/celda/decontx'
include { SCVITOOLS_SOLO               } from '../../modules/local/scvitools/solo'
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
    ch_h5ad = ch_datasets.h5ad.mix(ADATA_READRDS.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READRDS.out.versions)

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

    // Empty droplet detection


    // Doublet detection
    if (params.doublet_detection == 'solo') {
        SCVITOOLS_SOLO(ch_h5ad)
        ch_h5ad = SCVITOOLS_SOLO.out.h5ad
        ch_versions = SCVITOOLS_SOLO.out.versions
    }


    QC_FILTERED(ch_h5ad)
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    emit:
    h5ad          = ch_h5ad

    multiqc_files = ch_multiqc_files
    versions      = ch_versions                     // channel: [ versions.yml ]
}
