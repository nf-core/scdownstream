include { SCVITOOLS_SOLO                          } from '../../../modules/nf-core/scvitools/solo'
include { SCANPY_SCRUBLET                         } from '../../../modules/nf-core/scanpy/scrublet'
include { DOUBLETDETECTION                        } from '../../../modules/nf-core/doubletdetection'
include { SCDBLFINDER                             } from '../../../modules/nf-core/scdblfinder'
include { CUSTOM_DOUBLETREMOVAL as DOUBLETREMOVAL } from '../../../modules/local/custom/doubletremoval'

workflow DOUBLET_DETECTION {
    take:
    ch_h5ad // channel: [ meta, h5ad ]
    methods //   value: list of strings
    threshold //   value: integer
    doublet_removal //   value: boolean
    scvi_max_epochs //   value: integer

    main:
    ch_multiqc_files = channel.empty()
    ch_results = channel.empty()
    ch_obs = channel.empty()

    if (methods.size() == 0) {
        log.info("DOUBLET_DETECTION: Not performed since no methods selected.")
    }
    else {
        ch_batch_col = ch_h5ad.map { meta, _h5ad -> meta.batch_col }

        if (methods.contains('solo')) {
            SCVITOOLS_SOLO(
                ch_h5ad,
                ch_batch_col,
                scvi_max_epochs ?: [],
            )
            ch_results = ch_results.mix(SCVITOOLS_SOLO.out.predictions.map { meta, predictions -> [meta, 'solo', predictions] })
        }

        if (methods.contains('scrublet')) {
            SCANPY_SCRUBLET(
                ch_h5ad,
                ch_batch_col,
            )
            ch_results = ch_results.mix(SCANPY_SCRUBLET.out.predictions.map { meta, predictions -> [meta, 'scrublet', predictions] })
            ch_obs = ch_obs.mix(SCANPY_SCRUBLET.out.scores)
        }

        if (methods.contains('doubletdetection')) {
            DOUBLETDETECTION(
                ch_h5ad
            )
            ch_results = ch_results.mix(DOUBLETDETECTION.out.predictions.map { meta, predictions -> [meta, 'doubletdetection', predictions] })
            ch_obs = ch_obs.mix(DOUBLETDETECTION.out.scores)
        }

        if (methods.contains('scdblfinder')) {
            SCDBLFINDER(
                ch_h5ad
            )
            ch_results = ch_results.mix(SCDBLFINDER.out.predictions.map { meta, predictions -> [meta, 'scdblfinder', predictions] })
            ch_obs = ch_obs.mix(SCDBLFINDER.out.scores)
        }

        ch_results = ch_results
            .groupTuple()
            .map { meta, detector_methods, predictions ->
                def records = [detector_methods, predictions].transpose().sort { record -> record[0] }
                def method_order = records.collect { record -> record[0] }
                if (method_order != methods.sort(false)) {
                    error("Unexpected doublet detector results for '${meta.id}': expected ${methods}, received ${detector_methods}")
                }
                [meta, method_order, records.collect { record -> record[1] }]
            }

        DOUBLETREMOVAL(
            ch_h5ad.join(ch_results, failOnDuplicate: true, failOnMismatch: true),
            threshold,
            doublet_removal,
        )

        ch_h5ad = DOUBLETREMOVAL.out.h5ad
        ch_obs = ch_obs.mix(DOUBLETREMOVAL.out.obs)
        ch_multiqc_files = ch_multiqc_files.mix(DOUBLETREMOVAL.out.multiqc_files)
    }

    emit:
    h5ad          = ch_h5ad // channel: [ meta, h5ad ]
    obs           = ch_obs // channel: [ meta, csv/parquet/pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
}
