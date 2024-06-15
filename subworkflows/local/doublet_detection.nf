include { SCVITOOLS_SOLO  } from '../../modules/local/scvitools/solo'
include { SCANPY_SCRUBLET } from '../../modules/local/scanpy/scrublet'


workflow DOUBLET_DETECTION {
    take:
    ch_h5ad
    
    main:
    ch_versions = Channel.empty()

    methods = params.integration_methods.split(',').collect{it.trim().toLowerCase()}

    if (methods.contains('solo')) {
        SCVITOOLS_SOLO(ch_h5ad)
        ch_h5ad = SCVITOOLS_SOLO.out.h5ad
        ch_versions = SCVITOOLS_SOLO.out.versions
    }

    if (methods.contains('scrublet')) {
        SCANPY_SCRUBLET(ch_h5ad)
        ch_h5ad = SCANPY_SCRUBLET.out.h5ad
        ch_versions = SCANPY_SCRUBLET.out.versions
    }
    
    emit:
    h5ad = ch_h5ad

    versions = ch_versions
}