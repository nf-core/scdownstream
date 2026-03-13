include { SCVITOOLS_SOLO   } from '../../../modules/nf-core/scvitools/solo'
include { SCANPY_SCRUBLET  } from '../../../modules/nf-core/scanpy/scrublet'
include { DOUBLETDETECTION } from '../../../modules/nf-core/doubletdetection'
include { SCDS             } from '../../../modules/local/doublet_detection/scds'
include { SCDBLFINDER      } from '../../../modules/local/doublet_detection/scdblfinder'
include { DOUBLET_REMOVAL  } from '../../../modules/local/doublet_detection/doublet_removal'

workflow DOUBLET_DETECTION {
    take:
    ch_h5ad         // channel: [ meta, h5ad ]
    methods         //   value: list of strings
    threshold       //   value: integer
    scvi_max_epochs //   value: integer

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_predictions = channel.empty()

    if (methods.size() == 0) {
        log.info("DOUBLET_DETECTION: Not performed since no methods selected.")
    } else {
        ch_batch_col = ch_h5ad.map { meta, _h5ad -> meta.batch_col }
        ch_h5ad_doublet_rate = ch_h5ad.map { meta, h5ad -> [meta, h5ad, meta.doublet_rate] }

        if (methods.contains('scds')) {
            SCDS (
                ch_h5ad
            )
            ch_predictions = ch_predictions.mix(SCDS.out.predictions)
            ch_versions = ch_versions.mix(SCDS.out.versions)
        }

        if (methods.contains('solo')) {
            SCVITOOLS_SOLO (
                ch_h5ad,
                ch_batch_col,
                scvi_max_epochs ?: []
            )
            ch_predictions = ch_predictions.mix(SCVITOOLS_SOLO.out.predictions)
            ch_versions = ch_versions.mix(SCVITOOLS_SOLO.out.versions)
        }

        if (methods.contains('scrublet')) {
            SCANPY_SCRUBLET (
                ch_h5ad,
                ch_batch_col
            )
            ch_predictions = ch_predictions.mix(SCANPY_SCRUBLET.out.predictions)
            ch_versions = ch_versions.mix(SCANPY_SCRUBLET.out.versions)
        }

        if (methods.contains('doubletdetection')) {
            DOUBLETDETECTION (
                ch_h5ad
            )
            ch_predictions = ch_predictions.mix(DOUBLETDETECTION.out.predictions)
            ch_versions = ch_versions.mix(DOUBLETDETECTION.out.versions)
        }

        if (methods.contains('scdblfinder')) {
            SCDBLFINDER (
                ch_h5ad_doublet_rate
            )
            ch_predictions = ch_predictions.mix(SCDBLFINDER.out.predictions)
            ch_versions = ch_versions.mix(SCDBLFINDER.out.versions)
        }

        DOUBLET_REMOVAL (
            ch_h5ad.join(ch_predictions.groupTuple()),
            threshold,
        )

        ch_h5ad = DOUBLET_REMOVAL.out.h5ad
        ch_multiqc_files = ch_multiqc_files.mix(DOUBLET_REMOVAL.out.multiqc_files)
        ch_versions = ch_versions.mix(DOUBLET_REMOVAL.out.versions)
    }

    emit:
    h5ad          = ch_h5ad          // channel: [ meta, h5ad ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    versions      = ch_versions      // channel: [ versions.yml ]
}
