include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA as EMPTY_DROPLET_REMOVAL } from '../../nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { ANNDATA_GETSIZE as GET_UNFILTERED_SIZE                                     } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_FILTERED_SIZE                                       } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_THRESHOLDED_SIZE                                    } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_DEDOUBLETED_SIZE                                    } from '../../../modules/nf-core/anndata/getsize'
include { SCANPY_PLOTQC as QC_RAW                                                    } from '../../../modules/local/scanpy/plotqc'
include { AMBIENT_RNA_REMOVAL                                                        } from '../ambient_rna_removal'
include { SCANPY_FILTER                                                              } from '../../../modules/local/scanpy/filter'
include { DOUBLET_DETECTION                                                          } from '../doublet_detection'
include { SCANPY_PLOTQC as QC_FILTERED                                               } from '../../../modules/local/scanpy/plotqc'
include { CUSTOM_COLLECTSIZES as COLLECT_SIZES                                       } from '../../../modules/local/custom/collectsizes'

workflow QUALITY_CONTROL {
    take:
    ch_h5ad // channel: [ meta, filtered, unfiltered ]
    ambient_removal_method // value: string
    doublet_detection_methods // value: list of strings

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_sizes = Channel.empty()

    GET_UNFILTERED_SIZE(
        ch_h5ad.map { meta, filtered, unfiltered -> [meta, unfiltered ?: filtered] },
        "cells",
    )
    ch_versions = ch_versions.mix(GET_UNFILTERED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_UNFILTERED_SIZE.out.size.map { meta, size -> [meta.id, 'unfiltered', (size.text ?: "0").toInteger()] }
    )

    ch_h5ad = ch_h5ad.branch { meta, filtered, unfiltered ->
        complete: filtered
        return [meta, filtered, unfiltered]
        needs_filtering: unfiltered
        return [meta, filtered, unfiltered]
        problematic: true
        return [meta, filtered, unfiltered]
    }

    ch_complete = ch_h5ad.complete
    ch_needs_filtering = ch_h5ad.needs_filtering

    EMPTY_DROPLET_REMOVAL(ch_needs_filtering.map { meta, _filtered, unfiltered -> [meta, unfiltered] })
    ch_versions = ch_versions.mix(EMPTY_DROPLET_REMOVAL.out.versions)

    ch_complete = ch_complete.mix(
        ch_needs_filtering.join(EMPTY_DROPLET_REMOVAL.out.h5ad).map { meta, _empty, unfiltered, filtered -> [meta, filtered, unfiltered] }
    )

    GET_FILTERED_SIZE(
        ch_complete.map { meta, filtered, _unfiltered -> [meta, filtered] },
        "cells",
    )
    ch_versions = ch_versions.mix(GET_FILTERED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_FILTERED_SIZE.out.size.map { meta, size -> [meta.id, 'filtered', (size.text ?: "0").toInteger()] }
    )

    QC_RAW(ch_complete.map { meta, filtered, _unfiltered -> [meta, filtered] })
    ch_multiqc_files = ch_multiqc_files.mix(QC_RAW.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_RAW.out.versions)

    AMBIENT_RNA_REMOVAL(ch_complete, ambient_removal_method)
    ch_h5ad = AMBIENT_RNA_REMOVAL.out.h5ad
    ch_versions = ch_versions.mix(AMBIENT_RNA_REMOVAL.out.versions)

    ch_filtering = ch_h5ad.multiMap { meta, h5ad ->
        h5ad: [meta, h5ad]
        min_genes: meta.min_genes ?: 0
        min_cells: meta.min_cells ?: 0
        min_counts_gene: meta.min_counts_gene ?: 0
        min_counts_cell: meta.min_counts_cell ?: 0
        max_mito_percentage: meta.max_mito_percentage ?: 100
    }
    SCANPY_FILTER(ch_filtering.h5ad, ch_filtering.min_genes, ch_filtering.min_cells, ch_filtering.min_counts_gene, ch_filtering.min_counts_cell, ch_filtering.max_mito_percentage, params.dynamic_filtering)
    ch_h5ad = SCANPY_FILTER.out.h5ad
    ch_versions = ch_versions.mix(SCANPY_FILTER.out.versions)

    GET_THRESHOLDED_SIZE(ch_h5ad, "cells")
    ch_versions = ch_versions.mix(GET_THRESHOLDED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_THRESHOLDED_SIZE.out.size.map { meta, size -> [meta.id, 'thresholded', (size.text ?: "0").toInteger()] }
    )

    DOUBLET_DETECTION(ch_h5ad, doublet_detection_methods, params.doublet_detection_threshold)
    ch_h5ad = DOUBLET_DETECTION.out.h5ad
    ch_multiqc_files = ch_multiqc_files.mix(DOUBLET_DETECTION.out.multiqc_files)
    ch_versions = ch_versions.mix(DOUBLET_DETECTION.out.versions)

    GET_DEDOUBLETED_SIZE(ch_h5ad, "cells")
    ch_versions = ch_versions.mix(GET_DEDOUBLETED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_DEDOUBLETED_SIZE.out.size.map { meta, size -> [meta.id, 'dedoubleted', (size.text ?: "0").toInteger()] }
    )

    QC_FILTERED(ch_h5ad)
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    ch_sizes = ch_sizes
        .collectFile(
            seed: "sample\tstate\tsize",
            newLine: true,
            name: "size_list.tsv",
        ) { sample, state, size -> "${sample}\t${state}\t${size}" }
        .map { file -> [[id: 'sizes'], file] }

    COLLECT_SIZES(ch_sizes)
    ch_versions = ch_versions.mix(COLLECT_SIZES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_SIZES.out.multiqc_files)

    emit:
    h5ad          = ch_h5ad // channel: [ meta, h5ad ]
    multiqc_files = ch_multiqc_files // channel: [ json ]
    versions      = ch_versions // channel: [ versions.yml ]
}
