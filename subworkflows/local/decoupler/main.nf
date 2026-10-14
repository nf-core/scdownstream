include { DECOUPLER_NETWORK as NETWORK } from '../../../modules/local/decoupler/network'
include { DECOUPLER_ACTIVITY as ACTIVITY } from '../../../modules/local/decoupler/activity'
include { DECOUPLER_ENRICHMENT as ENRICHMENT } from '../../../modules/local/decoupler/enrichment'

workflow DECOUPLER {
    take:
    ch_pseudobulk_h5ad // channel: [ meta, anndata ]
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

    NETWORK(ch_network_requests, species)

    ch_networks = NETWORK.out.network
        .map { _meta, tsv -> tsv }
        .collect(sort: true)

    ACTIVITY(ch_pseudobulk_h5ad, ch_networks, method, tmin)

    ENRICHMENT(ch_de_results, ch_networks, method, tmin)

    emit:
    networks = NETWORK.out.network
    scores = ACTIVITY.out.scores
    enrichment = ENRICHMENT.out.results
    uns = ACTIVITY.out.uns.mix(ENRICHMENT.out.uns)
    multiqc_files = ACTIVITY.out.multiqc_files.mix(ENRICHMENT.out.multiqc_files).flatten()
}
