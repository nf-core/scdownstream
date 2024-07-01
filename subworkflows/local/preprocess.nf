include { SCANPY_READH5                } from '../../modules/local/scanpy/readh5'
include { ADATA_READRDS                } from '../../modules/local/adata/readrds'
include { ADATA_READCSV                } from '../../modules/local/adata/readcsv'
include { ADATA_UNIFY                  } from '../../modules/local/adata/unify'
include { SCANPY_PLOTQC as QC_RAW      } from '../../modules/local/scanpy/plotqc'
include { AMBIENT_RNA_REMOVAL          } from './ambient_rna_removal'
include { SCANPY_FILTER                } from '../../modules/local/scanpy/filter'
include { DOUBLET_DETECTION            } from './doublet_detection'
include { SCANPY_PLOTQC as QC_FILTERED } from '../../modules/local/scanpy/plotqc'

workflow PREPROCESS {

    take:
    ch_samples // channel: [ val(meta), file ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_h5ad = Channel.empty()

    ch_files = ch_samples.map { meta, file, raw -> [meta + [type: 'sample'], file] }
    ch_files = ch_files.mix(ch_samples
            .map { meta, file, raw -> [meta + [type: 'raw'], raw] }
            .filter { meta, raw -> raw }
        )

    ch_files = ch_files.map { meta, file -> [meta, file, file.extension.toLowerCase()] }
        .branch { meta, file, ext ->
            h5ad: ext == "h5ad"
                return [meta, file]
            h5: ext == "h5"
                return [meta, file]
            rds: ext == "rds"
                return [meta, file]
            csv: ext == "csv"
                return [meta, file]
        }

    ch_h5ad = ch_h5ad.mix(ch_files.h5ad)

    SCANPY_READH5(ch_files.h5)
    ch_h5ad = ch_h5ad.mix(SCANPY_READH5.out.h5ad)
    ch_versions = ch_versions.mix(SCANPY_READH5.out.versions)

    ADATA_READRDS(ch_files.rds)
    ch_h5ad = ch_h5ad.mix(ADATA_READRDS.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READRDS.out.versions)

    ADATA_READCSV(ch_files.csv)
    ch_h5ad = ch_h5ad.mix(ADATA_READCSV.out.h5ad)
    ch_versions = ch_versions.mix(ADATA_READCSV.out.versions)

    ch_samples = ch_h5ad.filter { meta, h5ad -> meta.type == 'sample' }
    ch_raws = ch_h5ad.filter { meta, h5ad -> meta.type == 'raw' }

    ADATA_UNIFY(ch_samples)
    ch_samples = ADATA_UNIFY.out.h5ad
    ch_versions = ch_versions.mix(ADATA_UNIFY.out.versions)

    QC_RAW(ch_samples)
    ch_multiqc_files = ch_multiqc_files.mix(QC_RAW.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_RAW.out.versions)

    AMBIENT_RNA_REMOVAL(ch_samples, ch_raws)
    ch_samples = AMBIENT_RNA_REMOVAL.out.h5ad
    ch_versions = ch_versions.mix(AMBIENT_RNA_REMOVAL.out.versions)

    SCANPY_FILTER(ch_samples)
    ch_samples = SCANPY_FILTER.out.h5ad
    ch_versions = ch_versions.mix(SCANPY_FILTER.out.versions)

    DOUBLET_DETECTION(ch_samples)
    ch_samples = DOUBLET_DETECTION.out.h5ad
    ch_multiqc_files = ch_multiqc_files.mix(DOUBLET_DETECTION.out.multiqc_files)
    ch_versions = ch_versions.mix(DOUBLET_DETECTION.out.versions)

    QC_FILTERED(ch_samples)
    ch_multiqc_files = ch_multiqc_files.mix(QC_FILTERED.out.multiqc_files)
    ch_versions = ch_versions.mix(QC_FILTERED.out.versions)

    emit:
    h5ad          = ch_samples

    multiqc_files = ch_multiqc_files
    versions      = ch_versions                     // channel: [ versions.yml ]
}
