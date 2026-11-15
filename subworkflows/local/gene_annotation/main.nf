include { ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION as GTF2FEATUREANNOTATION } from '../../../modules/nf-core/atlasgeneannotationmanipulation/gtf2featureannotation'
include { ADATA_GENEANNOTATION                                                          } from '../../../modules/local/adata/geneannotation'

workflow GENE_ANNOTATION {
    take:
    ch_h5ad // channel: [ meta, h5ad ]
    gtf     //    path: file or []
    species //   value: string

    main:
    ch_gene_table = channel.value([])

    if (gtf) {
        GTF2FEATUREANNOTATION(
            channel.value([[id: gtf.baseName], gtf]),
            channel.value([[], []]),
        )
        ch_gene_table = GTF2FEATUREANNOTATION.out.feature_annotation
            .map { _meta, tsv -> tsv }
            .first()
    }

    ADATA_GENEANNOTATION(
        ch_h5ad,
        ch_gene_table,
        species,
    )

    emit:
    var = ADATA_GENEANNOTATION.out.var.map { _meta, pkl -> pkl } // channel: [ pkl ]
    csv = ADATA_GENEANNOTATION.out.csv.map { _meta, csv -> csv } // channel: [ csv ]
}
