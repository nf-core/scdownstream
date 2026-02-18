#!/usr/bin/env python

import argparse
import sys
import scanpy as sc
import phate
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Run PHATE on an AnnData object.')
    parser.add_argument('--input', required=True, help='Input AnnData file (.h5ad)')
    parser.add_argument('--output', required=True, help='Output AnnData file (.h5ad)')
    parser.add_argument('--k', type=int, default=5, help='Number of nearest neighbors')
    parser.add_argument('--a', type=float, default=40, help='Alpha decay parameter')
    parser.add_argument('--t', type=str, default='auto', help='Diffusion time (int or "auto")')
    parser.add_argument('--n_jobs', type=int, default=1, help='Number of threads')
    parser.add_argument('--gamma', type=float, default=1, help='Informational distance parameter')
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Load data
    print(f"Reading input from {args.input}...")
    try:
        adata = sc.read_h5ad(args.input)
    except Exception as e:
        sys.exit(f"Error loading AnnData: {e}")

    # 2. Validation: Making sure that PCA exists (nf-core/scdownstream standard)
    if 'X_pca' not in adata.obsm:
        sys.exit("Error: 'X_pca' not found in adata.obsm. PHATE requires pre-computed PCA coordinates.")

    # 3. Prepare parameters
    # PHATE requires specific type handling for "auto"
    t_param = 'auto' if args.t == 'auto' else int(args.t)
    
    print(f"Running PHATE with k={args.k}, a={args.a}, t={t_param} on X_pca...")

    # 4. Run PHATE
    # We initialize the PHATE operator explicitly
    phate_op = phate.PHATE(
        n_pca=None,      # PCA is already computed
        knn=args.k,
        decay=args.a,
        t=t_param,
        gamma=args.gamma,
        n_jobs=args.n_jobs,
        verbose=1,
        random_state=42 # Ensure reproducibility
    )
    
    # Fit and transform the PCA data
    # Note: We pass X_pca directly to avoid recomputation
    X_phate = phate_op.fit_transform(adata.obsm['X_pca'])

    # 5. Store results
    # Save coordinates in the standard Scanpy slot
    adata.obsm['X_phate'] = X_phate

    # Save metadata for reproducibility and downstream usage (e.g. TDA)
    adata.uns['phate_params'] = {
        'k': args.k,
        'a': args.a,
        't': t_param, 
        'gamma': args.gamma,
        'diff_potential': phate_op.diff_potential # Crucial for Potential Distance
    }

    # 6. Save output
    print(f"Saving results to {args.output}...")
    try:
        adata.write_h5ad(args.output)
    except Exception as e:
        sys.exit(f"Error saving AnnData: {e}")

    print("PHATE computation completed successfully.")

if __name__ == "__main__":
    main()
