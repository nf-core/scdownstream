# EXPIMAP Pathway Database

## Overview

This directory contains pathway databases used by the EXPIMAP module for interpretable embedding and pathway analysis.

## Files

### pathways.gmt

- **Format**: Gene Matrix Transposed (GMT)
- **Content**: Reactome pathway gene sets
- **Description**: Contains biological pathway gene sets from the Reactome database that are used to construct the latent space in EXPIMAP models. Each row represents a pathway with its associated genes.
- **Source**: [Reactome Database](https://reactome.org/)
- **Version/Date**: Please refer to Reactome documentation for version information
- **Usage**: Used by the EXPIMAP module to guide the learning of interpretable latent representations in single-cell transcriptomics data

## Usage in Pipeline

The pathway file is automatically used by the EXPIMAP integration method. Users can specify a custom pathway file via the `expimap_gmt` parameter in the pipeline parameters, or use the default provided in this directory.

Example:

```bash
nextflow run nf-core/scdownstream \
  --input samplesheet.csv \
  --outdir results \
  --integration_methods expimap \
  --expimap_gmt assets/databases/expimap/pathways.gmt
```

## Notes

- The GMT format contains tab-separated values where each line represents one pathway
- First column: pathway name
- Second column: pathway URL/description
- Remaining columns: genes in the pathway
- The default Reactome pathways are suitable for most analyses of mammalian cells (especially human and mouse)
- For custom pathway analysis, users can provide their own GMT files using the `--expimap_gmt` parameter
