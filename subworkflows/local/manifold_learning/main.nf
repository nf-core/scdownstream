// subworkflows/local/manifold_learning/main.nf

include { PHATE     } from '../../../modules/local/phate/main'
include { DIFFMAP   } from '../../../modules/local/diffmap/main'
include { TOPOLOGY  } from '../../../modules/local/topology/main'
include { CAUSALITY } from '../../../modules/local/causality/main'

workflow MANIFOLD_LEARNING {

    take:
    ch_h5ad // Channel: [ val(meta), path(h5ad) ]
    methods // Value: String (e.g. "phate,diffmap")

    main:
    ch_versions = Channel.empty()
    ch_outputs  = Channel.empty()

    // ------------------------------------------------
    // 1. GEOMETRY STEP (PHATE / DIFFMAP)
    // ------------------------------------------------
    
    ch_geometry_out = Channel.empty()

    // Run PHATE if requested
    if ( methods.toString().toLowerCase().contains('phate') ) {
        PHATE ( ch_h5ad )
        ch_geometry_out = ch_geometry_out.mix( PHATE.out.h5ad )
        ch_versions     = ch_versions.mix( PHATE.out.versions )
    }

    // Run DIFFMAP if requested
    if ( methods.toString().toLowerCase().contains('diffmap') ) {
        DIFFMAP ( ch_h5ad )
        ch_geometry_out = ch_geometry_out.mix( DIFFMAP.out.h5ad )
        ch_versions     = ch_versions.mix( DIFFMAP.out.versions )
    }

    // If no geometry method ran (user error?), pass input through (fallback)
    // But ideally, we continue with the output of geometry
    
    // ------------------------------------------------
    // 2. TOPOLOGY & CAUSALITY STEPS
    // ------------------------------------------------
    // We run these on the output of the Geometry step.
    // Since PHATE/DIFFMAP output separate files, Topology/Causality will run for each.

    // Run Topology
    TOPOLOGY ( ch_geometry_out )
    ch_versions = ch_versions.mix( TOPOLOGY.out.versions )
    
    // Run Causality (Using Topology output to chain the information)
    CAUSALITY ( TOPOLOGY.out.h5ad )
    ch_versions = ch_versions.mix( CAUSALITY.out.versions )

    // Collect final outputs
    ch_outputs = CAUSALITY.out.h5ad

    emit:
    h5ad     = ch_outputs
    versions = ch_versions
}
