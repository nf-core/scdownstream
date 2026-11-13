include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA } from '../../nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { DROPLETUTILS_EMPTYDROPS                           } from '../../../modules/local/dropletutils/emptydrops'
include { ANNDATA_BARCODES                                  } from '../../../modules/nf-core/anndata/barcodes'

workflow EMPTY_DROPLET_REMOVAL {
    take:
    ch_unfiltered // channel: [ meta, h5ad ]
    method        //   value: string

    main:

    if (method == 'cellbender') {
        H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA (
            ch_unfiltered
        )
        ch_h5ad = H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA.out.h5ad
    }
    else if (method == 'emptydrops') {
        DROPLETUTILS_EMPTYDROPS (
            ch_unfiltered,
            params.emptydrops_lower,
            params.emptydrops_fdr
        )
        ANNDATA_BARCODES (
            ch_unfiltered.join(DROPLETUTILS_EMPTYDROPS.out.barcodes, failOnMismatch: true)
        )
        ch_h5ad = ANNDATA_BARCODES.out.h5ad
    }
    else {
        error("EMPTY_DROPLET_REMOVAL: Unexpected method for empty droplet removal: '${method}'.")
    }

    emit:
    h5ad = ch_h5ad // channel: [ meta, h5ad ]
}
