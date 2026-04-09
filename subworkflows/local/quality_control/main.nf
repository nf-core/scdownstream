include { SCANPY_CELLCYCLE                                                           } from '../../../modules/local/scanpy/cellcycle'
include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA as EMPTY_DROPLET_REMOVAL } from '../../nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { ANNDATA_GETSIZE as GET_UNFILTERED_SIZE                                     } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_FILTERED_SIZE                                       } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_THRESHOLDED_SIZE                                    } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_DEDOUBLETED_SIZE                                    } from '../../../modules/nf-core/anndata/getsize'
include { ANNDATA_GETSIZE as GET_SAMPLED_SIZE                                        } from '../../../modules/nf-core/anndata/getsize'
include { SCANPY_PLOTQC as QC_RAW                                                    } from '../../../modules/local/scanpy/plotqc'
include { AMBIENT_CORRECTION                                                         } from '../ambient_correction'
include { UNIFY                                                                      } from '../unify'
include { SCANPY_FILTER                                                              } from '../../../modules/local/scanpy/filter'
include { SCANPY_SAMPLE                                                              } from '../../../modules/local/scanpy/sample'
include { DOUBLET_DETECTION                                                          } from '../doublet_detection'
include { SCANPY_PLOTQC as QC_FILTERED                                               } from '../../../modules/local/scanpy/plotqc'
include { CUSTOM_COLLECTSIZES as COLLECT_SIZES                                       } from '../../../modules/local/custom/collectsizes'

workflow QUALITY_CONTROL {
    take:
    ch_h5ad                       // channel: [ meta, filtered, unfiltered ]
    ambient_correction_method     //   value: string
    ambient_corrected_integration //   value: boolean
    unify_gene_symbols            //   value: boolean
    duplicate_var_resolution      //   value: string
    aggregate_isoforms            //   value: boolean
    doublet_detection_methods     //   value: list of strings
    doublet_detection_threshold   //   value: float
    scvi_max_epochs               //   value: integer
    mito_genes                    //   value: string (path) or null
    sample_n                      //   value: string (integer > 1 or null)
    sample_fraction               //   value: string (float between 0-1 or null)
    cell_cycle_scoring            //   value: boolean
    s_genes                       //    path: file or []
    g2m_genes                     //    path: file or []

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_sizes = channel.empty()
    ch_obs_per_sample = channel.empty()

    GET_UNFILTERED_SIZE (
        ch_h5ad
        .map {
            meta, filtered, unfiltered ->
            [meta, unfiltered ?: filtered] },
        "cells",
    )
    ch_versions = ch_versions.mix(GET_UNFILTERED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_UNFILTERED_SIZE.out.size
        .map {
            meta, size ->
            [meta.id, 'unfiltered', (size.text ?: "0").toInteger()]
        }
    )

    ch_h5ad = ch_h5ad
        .branch {
            meta, filtered, unfiltered ->
            complete: filtered
            return [meta, filtered, unfiltered]
            needs_filtering: unfiltered
            return [meta, filtered, unfiltered]
            problematic: true
            return [meta, filtered, unfiltered]
        }

    ch_complete = ch_h5ad.complete
    ch_needs_filtering = ch_h5ad.needs_filtering

    EMPTY_DROPLET_REMOVAL (
        ch_needs_filtering
        .map {
            meta, _filtered, unfiltered ->
            [meta, unfiltered]
        }
    )

    ch_complete = ch_complete.mix(
        ch_needs_filtering
        .join(EMPTY_DROPLET_REMOVAL.out.h5ad)
        .map {
            meta, _empty, unfiltered, filtered ->
            [meta, filtered, unfiltered]
        }
    )

    GET_FILTERED_SIZE (
        ch_complete
            .map {
                meta, filtered, _unfiltered ->
                [meta, filtered]
            },
        "cells",
    )
    ch_versions = ch_versions.mix(GET_FILTERED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_FILTERED_SIZE.out.size
        .map {
            meta, size ->
            [meta.id, 'filtered', (size.text ?: "0").toInteger()]
        }
    )

    QC_RAW (
        ch_complete.map { meta, filtered, _unfiltered -> [meta, filtered] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(QC_RAW.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_RAW.out.versions)

    AMBIENT_CORRECTION (
        ch_complete,
        ambient_correction_method,
        ambient_corrected_integration
    )
    ch_h5ad = AMBIENT_CORRECTION.out.h5ad
    ch_versions = ch_versions.mix(AMBIENT_CORRECTION.out.versions)

    // Unification needds to happen before filtering to make sure all genes have symbols
    // Otherwise, mitochondrial gene detection will not work correctly
    UNIFY (
        ch_h5ad,
        unify_gene_symbols,
        duplicate_var_resolution,
        aggregate_isoforms
    )
    ch_versions = ch_versions.mix(UNIFY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(UNIFY.out.multiqc_files)
    ch_h5ad = UNIFY.out.h5ad

    ch_filtering = ch_h5ad
        .multiMap {
            meta, h5ad ->
            h5ad: [meta, h5ad]
            symbol_col: meta.symbol_col ?: "index"
            min_genes: meta.min_genes ?: 0
            min_cells: meta.min_cells ?: 0
            min_counts_gene: meta.min_counts_gene ?: 0
            min_counts_cell: meta.min_counts_cell ?: 0
            max_mito_percentage: meta.max_mito_percentage ?: 100
        }
    SCANPY_FILTER (
        ch_filtering.h5ad,
        ch_filtering.symbol_col,
        ch_filtering.min_genes,
        ch_filtering.min_cells,
        ch_filtering.min_counts_gene,
        ch_filtering.min_counts_cell,
        ch_filtering.max_mito_percentage,
        mito_genes ?: []
    )
    ch_h5ad = SCANPY_FILTER.out.h5ad
    ch_versions = ch_versions.mix(SCANPY_FILTER.out.versions)

    // Only run SCANPY_SAMPLE if sample_n or sample_fraction is set
    if (sample_n || sample_fraction) {
        SCANPY_SAMPLE (
            ch_h5ad,
            sample_n ?: [],
            sample_fraction ?: []
        )
        ch_h5ad = SCANPY_SAMPLE.out.h5ad
        ch_versions = ch_versions.mix(SCANPY_SAMPLE.out.versions)

        GET_SAMPLED_SIZE (
            ch_h5ad,
            "cells"
        )
        ch_versions = ch_versions.mix(GET_SAMPLED_SIZE.out.versions)
        ch_sizes = ch_sizes.mix(
            GET_SAMPLED_SIZE.out.size
            .map {
                meta, size ->
                [meta.id, 'sampled', (size.text ?: "0").toInteger()]
            }
        )
    }

    GET_THRESHOLDED_SIZE (
        ch_h5ad,
        "cells"
    )
    ch_versions = ch_versions.mix(GET_THRESHOLDED_SIZE.out.versions)
    ch_sizes = ch_sizes.mix(
        GET_THRESHOLDED_SIZE.out.size
        .map {
            meta, size ->
            [meta.id, 'thresholded', (size.text ?: "0").toInteger()]
        }
    )

    DOUBLET_DETECTION (
        ch_h5ad,
        doublet_detection_methods,
        doublet_detection_threshold,
        scvi_max_epochs
    )
    ch_h5ad = DOUBLET_DETECTION.out.h5ad
    ch_multiqc_files = ch_multiqc_files.mix(DOUBLET_DETECTION.out.multiqc_files)
    ch_versions = ch_versions.mix(DOUBLET_DETECTION.out.versions)

    if (doublet_detection_methods.size() > 0) {
        GET_DEDOUBLETED_SIZE (
            ch_h5ad,
            "cells"
        )
        ch_versions = ch_versions.mix(GET_DEDOUBLETED_SIZE.out.versions)
        ch_sizes = ch_sizes.mix(
            GET_DEDOUBLETED_SIZE.out.size
            .map {
                meta, size ->
                [meta.id, 'dedoubleted', (size.text ?: "0").toInteger()]
            }
        )
    }

    QC_FILTERED (
        ch_h5ad
    )
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    if (cell_cycle_scoring) {
        ch_cellcycle = ch_h5ad.multiMap {
            meta, h5ad ->
            h5ad:       [meta, h5ad]
            symbol_col: meta.symbol_col ?: "index"
        }
        SCANPY_CELLCYCLE (
            ch_cellcycle.h5ad,
            s_genes,
            g2m_genes,
            ch_cellcycle.symbol_col
        )
        ch_obs_per_sample = ch_obs_per_sample.mix(SCANPY_CELLCYCLE.out.obs)
        ch_versions = ch_versions.mix(SCANPY_CELLCYCLE.out.versions)
    }

    ch_sizes = ch_sizes
        .collectFile(
            seed: "sample\tstate\tsize",
            newLine: true,
            name: "size_list.tsv",
        ) { sample, state, size -> "${sample}\t${state}\t${size}" }
        .map { file -> [[id: 'sizes'], file] }

    COLLECT_SIZES (
        ch_sizes
    )
    ch_versions = ch_versions.mix(COLLECT_SIZES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_SIZES.out.multiqc_files)

    emit:
    h5ad          = ch_h5ad           // channel: [ meta, h5ad ]
    sizes         = ch_sizes          // channel: [ tsv ]
    obs           = ch_obs_per_sample // channel: [ meta, pkl ]
    multiqc_files = ch_multiqc_files  // channel: [ json ]
    versions      = ch_versions       // channel: [ versions.yml ]
}
