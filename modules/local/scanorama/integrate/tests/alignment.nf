include { SCANORAMA_INTEGRATE                       } from '../main'
include { ADATA_EXTEND as EXTEND_SCANORAMA_FRAGMENT } from '../../../adata/extend/main'

workflow SCANORAMA_ALIGNMENT {
    take:
    h5ad

    main:
    SCANORAMA_INTEGRATE(h5ad, 'regression_batch', 'X', true)
    fragments = SCANORAMA_INTEGRATE.out.obsm.map { fragment -> [fragment.baseName.substring(2), fragment] }
    EXTEND_SCANORAMA_FRAGMENT(
        h5ad.map { meta, base -> [meta.id, meta, base] }.join(fragments, failOnDuplicate: true, failOnMismatch: true).map { _id, meta, base, fragment -> [meta, base, [], [], [fragment], [], [], []] }
    )

    emit:
    integrated = SCANORAMA_INTEGRATE.out.h5ad
    extended   = EXTEND_SCANORAMA_FRAGMENT.out.h5ad
    versions   = SCANORAMA_INTEGRATE.out.versions
}
