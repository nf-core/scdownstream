include { DECOUPLER_NETWORK    } from '../../../modules/local/decoupler/network'
include { DECOUPLER_ACTIVITY   } from '../../../modules/local/decoupler/activity'
include { DECOUPLER_ENRICHMENT } from '../../../modules/local/decoupler/enrichment'

workflow PATHWAY_ACTIVITY {
    take:
    ch_pseudobulk_h5ad // channel: [ meta, h5ad ]
    ch_de_results // channel: [ meta, [ parquet ] ] with meta.de_method
    resources //   value: string, comma-separated Omnipath resources
    custom_network //    path: custom TSV/GMT network or []
    species //   value: string
    method //   value: string
    tmin //   value: integer

    main:
    def network_requests = (resources ?: '').tokenize(',')*.trim().findAll { name -> name }.unique().collect { name -> [[id: name], []] }
    if (custom_network) {
        network_requests << [[id: custom_network.baseName], custom_network]
    }

    // Avoids Omnipath queries when no grouping requests decoupler
    ch_network_requests = ch_pseudobulk_h5ad
        .mix(ch_de_results)
        .first()
        .flatMap { _item -> network_requests }

    DECOUPLER_NETWORK(ch_network_requests, species)

    ch_networks = DECOUPLER_NETWORK.out.network
        .map { _meta, tsv -> tsv }
        .collect(sort: true)

    DECOUPLER_ACTIVITY(ch_pseudobulk_h5ad, ch_networks, method, tmin)

    DECOUPLER_ENRICHMENT(ch_de_results, ch_networks, method, tmin)

    emit:
    networks      = DECOUPLER_NETWORK.out.network
    scores        = DECOUPLER_ACTIVITY.out.scores
    enrichment    = DECOUPLER_ENRICHMENT.out.results
    uns           = DECOUPLER_ACTIVITY.out.uns.mix(DECOUPLER_ENRICHMENT.out.uns)
    multiqc_files = DECOUPLER_ACTIVITY.out.multiqc_files.mix(DECOUPLER_ENRICHMENT.out.multiqc_files).flatten()
}
