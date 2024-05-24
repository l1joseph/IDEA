import numpy as np
import anndata as ad

def median_of_ratios(adata):
    """
    Calculate size factors using the median-of-ratios method.
    
    Parameters:
    adata (anndata.AnnData): AnnData object containing the raw counts matrix.

    Returns:
    np.ndarray: Array of size factors for each sample.
    """
    # Calculate the geometric mean of each gene across all samples
    geometric_means = np.exp(np.log(adata.X + 1).mean(axis=0))
    
    # Compute ratios for each gene in each sample to the geometric mean
    ratios = adata.X / geometric_means
    
    # Median of ratios for each sample
    size_factors = np.median(ratios, axis=1)
    
    return size_factors

def size_factors(adata):
    """
    Calculate and store size factors in the AnnData object.
    
    Parameters:
    adata (anndata.AnnData): AnnData object containing the raw counts matrix.

    Returns:
    pd.Series: Pandas Series of size factors with sample identifiers.
    """
    adata.obs['size_factors'] = median_of_ratios(adata)
    return adata.obs['size_factors']

def normalize_counts(adata):
    """
    Normalize the counts matrix in the AnnData object using size factors.
    
    Parameters:
    adata (anndata.AnnData): AnnData object containing the raw counts matrix.

    Returns:
    anndata.AnnData: AnnData object with normalized counts matrix stored in `adata.layers['normalized']`.
    """
    if 'size_factors' not in adata.obs:
        size_factors(adata)
    
    # Normalize counts
    adata.layers['normalized'] = adata.X / adata.obs['size_factors'].values[:, np.newaxis]
    return adata
