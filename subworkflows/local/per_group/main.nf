include { SCANPY_PAGA            } from '../../../modules/local/scanpy/paga'
include { SCANPY_RANKGENESGROUPS } from '../../../modules/local/scanpy/rankgenesgroups'
include { LIANA_RANKAGGREGATE    } from '../../../modules/local/liana/rankaggregate'
include { DIFFERENTIAL_EXPRESSION } from '../differential_expression'
include { CYTETYPE                 } from '../../../modules/local/cytetype'

workflow PER_GROUP {
    take:
    ch_h5ad_with_neighbors // channel: [ meta, h5ad ], anndata objects with neighbors, one per embedding and annotation
    ch_h5ad_no_neighbors   // channel: [ meta, h5ad ], anndata objects without neighbors, one per annotation
    skip_liana             //   value: boolean
    skip_rankgenesgroups   //   value: boolean
    cytetype_study_context //   value: string

    main:
    ch_uns           = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_h5ad          = channel.empty()
    ch_obs           = channel.empty()

    SCANPY_PAGA(
        ch_h5ad_with_neighbors
            .filter { meta, _h5ad -> meta.analyses == null || 'paga' in meta.analyses }
    )
    ch_uns           = ch_uns.mix(SCANPY_PAGA.out.uns)
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_PAGA.out.multiqc_files)

    if (!skip_liana) {
        LIANA_RANKAGGREGATE(
            ch_h5ad_no_neighbors
                .filter { meta, _h5ad -> meta.analyses == null || 'liana' in meta.analyses }
        )
        ch_uns           = ch_uns.mix(LIANA_RANKAGGREGATE.out.uns)
    }

    if (!skip_rankgenesgroups) {
        ch_h5ad_for_de = ch_h5ad_no_neighbors
            .filter { meta, _h5ad -> meta.analyses == null || 'de' in meta.analyses }

        DIFFERENTIAL_EXPRESSION(
            ch_h5ad_for_de
        )
        ch_uns           = ch_uns.mix(DIFFERENTIAL_EXPRESSION.out.uns)
        ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_EXPRESSION.out.multiqc_files)
        ch_h5ad          = ch_h5ad.mix(DIFFERENTIAL_EXPRESSION.out.h5ad)

        if (cytetype_study_context) {
            ch_h5ad_for_cytetype = ch_h5ad
                .filter { meta, _h5ad -> meta.analyses == null || 'cytetype' in meta.analyses }

            CYTETYPE(
                ch_h5ad_for_cytetype,
                "index",
                cytetype_study_context,
                ch_h5ad_for_cytetype.map { meta, _h5ad -> meta.obs_key },
                "rank_genes_groups"
            )
            ch_obs = ch_obs.mix(CYTETYPE.out.obs)
        }
    }

    emit:
    uns           = ch_uns           // channel: [ pkl ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    h5ad          = ch_h5ad            // channel: [ meta, h5ad ] — empty if skip_rankgenesgroups
    obs           = ch_obs             // channel: [ meta, pkl ] — empty if cytetype skipped
}
