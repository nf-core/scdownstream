include { SCVITOOLS_SOLO  } from '../../modules/local/scvitools/solo'
include { SCANPY_SCRUBLET } from '../../modules/local/scanpy/scrublet'
include { DOUBLET_REMOVAL } from '../../modules/local/doublet_removal'

workflow DOUBLET_DETECTION {
    take:
    ch_h5ad
    
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_predictions = Channel.empty()

    methods = params.doublet_detection.split(',').collect{it.trim().toLowerCase()}

    if (methods.contains('solo')) {
        SCVITOOLS_SOLO(ch_h5ad.map{meta, h5ad -> [[id: 'solo'], h5ad]})
        ch_predictions = ch_predictions.mix(SCVITOOLS_SOLO.out.predictions)
        ch_versions = SCVITOOLS_SOLO.out.versions
    }

    if (methods.contains('scrublet')) {
        SCANPY_SCRUBLET(ch_h5ad.map{meta, h5ad -> [[id: 'scrublet'], h5ad]})
        ch_predictions = ch_predictions.mix(SCANPY_SCRUBLET.out.predictions)
        ch_versions = SCANPY_SCRUBLET.out.versions
    }

    DOUBLET_REMOVAL(
        ch_h5ad,
        ch_predictions.flatten().collect(),
        params.doublet_detection_threshold
    )

    ch_h5ad = DOUBLET_REMOVAL.out.h5ad
    ch_versions = ch_versions.mix(DOUBLET_REMOVAL.out.versions)

    emit:
    h5ad = ch_h5ad

    multiqc_files = ch_multiqc_files
    versions = ch_versions
}
