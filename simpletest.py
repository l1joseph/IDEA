import pandas as pd
import anndata as ad
from dataset import DESeqDataSet
from normalization import size_factors, normalize_counts

def main():
    # Sample data
    counts = pd.DataFrame({
        'gene1': [10, 15, 20],
        'gene2': [5, 8, 12],
        'gene3': [30, 45, 50]
    }, index=['sample1', 'sample2', 'sample3'])

    sample_info = pd.DataFrame({
        'condition': ['A', 'B', 'A']
    }, index=['sample1', 'sample2', 'sample3'])

    # Convert to AnnData
    adata = ad.AnnData(X=counts.values, obs=sample_info, var=pd.DataFrame(index=counts.columns))

    # Print initial data
    print("Initial Counts:\n", adata.to_df())
    print("Sample Info:\n", adata.obs)

    # Perform normalization
    size_factors(adata)
    normalize_counts(adata)

    # Print results
    print("Size Factors:\n", adata.obs['size_factors'])
    print("Normalized Counts:\n", pd.DataFrame(adata.layers['normalized'], index=adata.obs.index, columns=adata.var.index))

if __name__ == "__main__":
    main()
