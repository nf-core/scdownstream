include { samplesheetToList    } from 'plugin/nf-schema'
include { SINGLER              } from '../singler'
include { CELLTYPES_CELLTYPIST } from '../../../modules/local/celltypist'

workflow PER_CELL_ANNOTATION {
    take:
    ch_h5ad                   // channel: [ meta, h5ad, symbol_col, counts_layer ]
    celldex_reference         //   value: string
    celltypist_model          //   value: string

    main:
    ch_obs = channel.empty()
    ch_annotation_column_rows = channel.empty()

    if (celldex_reference ) {
        SINGLER (
            ch_h5ad,
            channel.fromList(samplesheetToList(
                celldex_reference,
                "${projectDir}/assets/schema_singler.json")
            )
        )
        ch_obs = ch_obs.mix(SINGLER.out.obs)
        ch_annotation_column_rows = ch_annotation_column_rows.mix(SINGLER.out.annotation_columns)
    }

    if (celltypist_model) {
        celltypist_models = channel.value(celltypist_model
            .split(',')
            .collect{ it -> it.trim() }
        )

        CELLTYPES_CELLTYPIST (
            ch_h5ad.map { meta, h5ad, symbol_col, _counts_layer -> [meta, h5ad, symbol_col] },
            celltypist_models
        )
        ch_obs = ch_obs.mix(CELLTYPES_CELLTYPIST.out.obs)
        ch_annotation_column_rows = ch_annotation_column_rows.mix(CELLTYPES_CELLTYPIST.out.annotation_columns)
    }

    ch_annotation_column_rows = ch_annotation_column_rows
        .splitCsv(header: true, elem: 1)

    emit:
    obs                    = ch_obs                    // channel: [ meta, parquet ]
    annotation_column_rows = ch_annotation_column_rows // channel: [ obs_column, aggregatable ]
}
