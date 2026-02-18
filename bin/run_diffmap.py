#!/usr/bin/env python

import argparse
import sys
import scanpy as sc
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Run Diffusion Maps (Diffmap) on an AnnData object.')
    parser.add_argument('--input', required=True, help='Input AnnData file (.h5ad)')
    parser.add_argument('--output', required=True, help='Output AnnData file (.h5ad)')
    parser.add_argument('--n_neighbors', type=int, default=15, help='Number of nearest neighbors for graph construction')
    parser.add_argument('--n_comps', type=int, default=15, help='Number of diffusion components to compute')
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Load data
    print(f"Reading input from {args.input}...")
    try:
        adata = sc.read_h5ad(args.input)
    except Exception as e:
        sys.exit(f"Error loading AnnData: {e}")

    # 2. Validation
    if 'X_pca' not in adata.obsm:
        sys.exit("Error: 'X_pca' not found in adata.obsm. Diffmap requires pre-computed PCA coordinates.")

    # 3. Compute Neighbors
    # Diffmap requires a neighborhood graph. We compute it on X_pca.
    print(f"Computing neighbors graph (k={args.n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep='X_pca')

    # 4. Run diffusion maps
    print(f"Running Diffusion Maps with n_comps={args.n_comps}...")
    try:
        sc.tl.diffmap(adata, n_comps=args.n_comps)
    except Exception as e:
        sys.exit(f"Error running diffmap: {e}")

    # The result is stored in adata.obsm['X_diffmap'] automatically by scanpy

    # 5. Save output
    print(f"Saving results to {args.output}...")
    try:
        adata.write_h5ad(args.output)
    except Exception as e:
        sys.exit(f"Error saving AnnData: {e}")

    # 6. Generate versions.yml
    import scanpy as s_lib
    
    with open("versions.yml", "w") as f:
        f.write('process:\n')
        f.write(f'  scanpy: "{s_lib.__version__}"\n')
        f.write(f'  numpy: "{np.__version__}"\n')

    print("Diffmap computation completed successfully.")

if __name__ == "__main__":
    main()
