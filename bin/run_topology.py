#!/usr/bin/env python

import argparse
import sys
import scanpy as sc
import numpy as np

# Secure import for ripser
try:
    from ripser import ripser
except ImportError:
    # Fail gracefully if the library is missing
    sys.exit("Error: 'ripser' library not found. Please ensure it is installed in the environment.")

def parse_args():
    parser = argparse.ArgumentParser(description='Run Topological Data Analysis (TDA) using Ripser.')
    parser.add_argument('--input', required=True, help='Input AnnData file (.h5ad)')
    parser.add_argument('--output', required=True, help='Output AnnData file (.h5ad)')
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Load Data
    print(f"Reading input from {args.input}...")
    try:
        adata = sc.read_h5ad(args.input)
    except Exception as e:
        sys.exit(f"Error loading AnnData: {e}")

    # 2. Select Embedding
    # Priority: PHATE -> Diffmap -> PCA
    # We check which embeddings are available in the object
    embedding_key = 'X_phate' 
    if embedding_key not in adata.obsm:
        if 'X_diffmap' in adata.obsm:
            embedding_key = 'X_diffmap'
        elif 'X_pca' in adata.obsm:
            embedding_key = 'X_pca'
        else:
            sys.exit("Error: No valid embedding found (requires X_phate, X_diffmap, or X_pca).")

    print(f"Running Ripser TDA on {embedding_key}...")
    
    # 3. Subsampling cause TDA is computationally expensive
    # If the dataset is large, we subsample to make sure the pipeline doesn't hang.
    matrix = adata.obsm[embedding_key]
    if matrix.shape[0] > 1000:
        print("Subsampling to 1000 cells for performance...")
        # Use random sampling without replacement
        idx = np.random.choice(matrix.shape[0], 1000, replace=False)
        matrix = matrix[idx, :]

    # 4. Run Persistent Homology
    try:
        # maxdim=1 computes H0 (connected components) and H1 (loops)
        diagrams = ripser(matrix, maxdim=1)['dgms']
    except Exception as e:
        sys.exit(f"Error running ripser: {e}")

    # 5. Store results
    diagrams_dict = {}
    for i, dgm in enumerate(diagrams):
        diagrams_dict[f'dim_{i}'] = dgm

    adata.uns['tda_results'] = {
        'max_homology_dim': 1, 
        'embedding_used': embedding_key,
        'diagrams': diagrams_dict 
    }

    # 6. Save output
    print(f"Saving results to {args.output}...")
    try:
        adata.write_h5ad(args.output)
    except Exception as e:
        sys.exit(f"Error saving h5ad file: {e}")
    
    # 7. Generate versions
    import ripser as r
    with open("versions.yml", "w") as f:
        f.write(f'process:\n')
        f.write(f'  ripser: "{r.__version__}"\n')
        f.write(f'  scanpy: "{sc.__version__}"\n')
        f.write(f'  numpy: "{np.__version__}"\n')

    print("TDA computation completed successfully.")

if __name__ == "__main__":
    main()
