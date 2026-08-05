include { SCANPY_CELLCYCLE                                                           } from '../../../modules/local/scanpy/cellcycle'
include { H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA as EMPTY_DROPLET_REMOVAL } from '../../nf-core/h5ad_removebackground_barcodes_cellbender_anndata'
include { SCANPY_PLOTQC as QC_RAW                                                    } from '../../../modules/local/scanpy/plotqc'
include { AMBIENT_CORRECTION                                                         } from '../ambient_correction'
include { UNIFY                                                                      } from '../unify'
include { SCANPY_FILTER                                                              } from '../../../modules/local/scanpy/filter'
include { SCANPY_SAMPLE                                                              } from '../../../modules/local/scanpy/sample'
include { DOUBLET_DETECTION                                                          } from '../doublet_detection'
include { SCANPY_PLOTQC as QC_FILTERED                                               } from '../../../modules/local/scanpy/plotqc'
include { CUSTOM_COLLECTSIZES as COLLECT_SIZES                                       } from '../../../modules/local/custom/collectsizes'
include { anndata                                                                    } from 'plugin/nf-anndata'

def countCells(h5ad) {
    workflow.stubRun ? 0 : anndata(h5ad).n_obs
}

workflow QUALITY_CONTROL {
    take:
    ch_h5ad                       // channel: [ meta, filtered, unfiltered ]
    ambient_correction_method     //   value: string
    ambient_corrected_integration //   value: boolean
    unify_gene_symbols            //   value: boolean
    duplicate_var_resolution      //   value: string
    aggregate_isoforms            //   value: boolean
    doublet_detection_methods     //   value: list of strings
    doublet_detection_threshold   //   value: integer
    doublet_removal               //   value: boolean
    scvi_max_epochs               //   value: integer
    mito_genes                    //   value: string (path) or null
    sample_n                      //   value: string (integer > 1 or null)
    sample_fraction               //   value: string (float between 0-1 or null)
    cell_cycle_scoring            //   value: boolean
    s_genes                       //    path: file or []
    g2m_genes                     //    path: file or []
    species                       //   value: string

    main:
    ch_multiqc_files = channel.empty()
    ch_sizes = channel.empty()
    ch_obs_per_sample = channel.empty()

    ch_sizes = ch_sizes.mix(
        ch_h5ad.map { meta, filtered, unfiltered ->
            [meta.id, 'unfiltered', countCells(unfiltered ?: filtered)]
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

    ch_sizes = ch_sizes.mix(
        ch_complete.map { meta, filtered, _unfiltered ->
            [meta.id, 'filtered', countCells(filtered)]
        }
    )

    ch_qc_plot = ch_complete.multiMap {
        meta, filtered, _unfiltered ->
        h5ad: [meta, filtered]
        symbol_col: meta.symbol_col ?: "index"
    }
    QC_RAW (
        ch_qc_plot.h5ad,
        ch_qc_plot.symbol_col,
        mito_genes ?: [],
        'Unfiltered QC plots',
        'Quality control plots',
    )
    ch_multiqc_files = ch_multiqc_files.mix(QC_RAW.out.multiqc_files)

    AMBIENT_CORRECTION (
        ch_complete,
        ambient_correction_method,
        ambient_corrected_integration
    )
    ch_h5ad = AMBIENT_CORRECTION.out.h5ad

    // Unification needds to happen before filtering to make sure all genes have symbols
    // Otherwise, mitochondrial gene detection will not work correctly
    UNIFY (
        ch_h5ad,
        unify_gene_symbols,
        duplicate_var_resolution,
        aggregate_isoforms,
        species
    )
    ch_multiqc_files = ch_multiqc_files.mix(UNIFY.out.multiqc_files)
    ch_h5ad = UNIFY.out.h5ad

    ch_filtering = ch_h5ad
        .multiMap {
            meta, h5ad ->
            h5ad: [meta, h5ad]
            symbol_col: meta.symbol_col ?: "index"
            min_genes: meta.min_genes
            min_cells: meta.min_cells
            min_counts_gene: meta.min_counts_gene
            min_counts_cell: meta.min_counts_cell
            max_mito_percentage: meta.max_mito_percentage
            min_ribo_percentage: meta.min_ribo_percentage
            max_hb_percentage: meta.max_hb_percentage
            log1p_total_counts_nmads: meta.log1p_total_counts_nmads
            log1p_n_genes_by_counts_nmads: meta.log1p_n_genes_by_counts_nmads
            pct_counts_in_top_20_genes_nmads: meta.pct_counts_in_top_20_genes_nmads
            pct_counts_mt_nmads: meta.pct_counts_mt_nmads
        }
    SCANPY_FILTER (
        ch_filtering.h5ad,
        ch_filtering.symbol_col,
        ch_filtering.min_genes,
        ch_filtering.min_cells,
        ch_filtering.min_counts_gene,
        ch_filtering.min_counts_cell,
        ch_filtering.max_mito_percentage,
        ch_filtering.min_ribo_percentage,
        ch_filtering.max_hb_percentage,
        ch_filtering.log1p_total_counts_nmads,
        ch_filtering.log1p_n_genes_by_counts_nmads,
        ch_filtering.pct_counts_in_top_20_genes_nmads,
        ch_filtering.pct_counts_mt_nmads,
        mito_genes ?: [],
        true,
        'Filter threshold histograms',
        'QC metric histograms with applied filter thresholds',
    )
    ch_h5ad = SCANPY_FILTER.out.h5ad
        .map { meta, h5ad ->
            if (!workflow.stubRun) {
                def ad = anndata(h5ad)
                if (ad.n_obs == 0) {
                    error("No cells remaining after filtering for sample '${meta.id}'")
                }
                if (ad.n_vars == 0) {
                    error("No genes remaining after filtering for sample '${meta.id}'")
                }
            }
            [meta, h5ad]
        }
    ch_multiqc_files = ch_multiqc_files.mix(SCANPY_FILTER.out.multiqc_files.flatten())

    // Only run SCANPY_SAMPLE if sample_n or sample_fraction is set
    if (sample_n || sample_fraction) {
        SCANPY_SAMPLE (
            ch_h5ad,
            sample_n ?: [],
            sample_fraction ?: []
        )
        ch_h5ad = SCANPY_SAMPLE.out.h5ad

        ch_sizes = ch_sizes.mix(
            ch_h5ad.map { meta, h5ad ->
                [meta.id, 'sampled', countCells(h5ad)]
            }
        )
    }

    ch_sizes = ch_sizes.mix(
        ch_h5ad.map { meta, h5ad ->
            [meta.id, 'thresholded', countCells(h5ad)]
        }
    )

    DOUBLET_DETECTION (
        ch_h5ad,
        doublet_detection_methods,
        doublet_detection_threshold,
        doublet_removal,
        scvi_max_epochs
    )
    ch_h5ad = DOUBLET_DETECTION.out.h5ad
    ch_multiqc_files = ch_multiqc_files.mix(DOUBLET_DETECTION.out.multiqc_files)

    if (doublet_detection_methods.size() > 0 && doublet_removal) {
        ch_sizes = ch_sizes.mix(
            ch_h5ad.map { meta, h5ad ->
                [meta.id, 'dedoubleted', countCells(h5ad)]
            }
        )
    }

    ch_qc_filtered_plot = ch_h5ad.multiMap {
        meta, h5ad ->
        h5ad: [meta, h5ad]
        symbol_col: meta.symbol_col ?: "index"
    }
    QC_FILTERED (
        ch_qc_filtered_plot.h5ad,
        ch_qc_filtered_plot.symbol_col,
        mito_genes ?: [],
        'Filtered QC plots',
        'Quality control plots',
    )
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)

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
        ch_multiqc_files = ch_multiqc_files.mix(SCANPY_CELLCYCLE.out.multiqc_files)
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
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_SIZES.out.multiqc_files)

    emit:
    h5ad          = ch_h5ad           // channel: [ meta, h5ad ]
    sizes         = ch_sizes          // channel: [ tsv ]
    obs           = ch_obs_per_sample // channel: [ meta, parquet ]
    multiqc_files = ch_multiqc_files  // channel: [ json ]
}
