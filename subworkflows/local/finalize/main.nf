include { ADATA_EXTEND } from '../../../modules/local/adata/extend'
include { ADATA_TORDS } from '../../../modules/local/adata/tords'
include { ADATA_PREPCELLXGENE } from '../../../modules/local/adata/prepcellxgene'

workflow FINALIZE {
    take:
    ch_h5ad // channel: [ merged, anndata ]
    ch_obs // channel: [ pkl ]
    ch_var // channel: [ pkl ]
    ch_obsm // channel: [ pkl ]
    ch_obsp
    ch_uns // channel: [ pkl ]
    prep_cellxgene //   value: boolean
    tords //   value: boolean
    skip_qc_report //   value: boolean

    main:

    ADATA_EXTEND(
        ch_h5ad.combine(ch_obs.flatten().collect().ifEmpty([]).map { it -> [it] }).combine(ch_var.flatten().collect().ifEmpty([]).map { it -> [it] }).combine(ch_obsm.flatten().collect().ifEmpty([]).map { it -> [it] }).combine(ch_obsp.flatten().collect().ifEmpty([]).map { it -> [it] }).combine(ch_uns.flatten().collect().ifEmpty([]).map { it -> [it] }).combine(channel.value([[]])),
        (!skip_qc_report || tords) ? 'h5ad' : 'zarr',
    )

    if (tords) {
        ADATA_TORDS(
            ADATA_EXTEND.out.anndata,
            'X',
        )
    }

    ch_h5ad_cellxgene = channel.empty()
    if (prep_cellxgene) {
        ADATA_PREPCELLXGENE(
            ADATA_EXTEND.out.anndata,
            'h5ad',
        )
        ch_h5ad_cellxgene = ADATA_PREPCELLXGENE.out.anndata
    }

    emit:
    anndata = ADATA_EXTEND.out.anndata // channel: [ meta, anndata ]
    anndata_cellxgene = ch_h5ad_cellxgene // channel: [ meta, h5ad ]
}
