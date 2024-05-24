import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def estimate_dispersions(dds):
    """
    Estimate the dispersions for each gene.
    
    Parameters:
    dds (anndata.AnnData): AnnData object containing the counts matrix and necessary metadata.

    Returns:
    pd.Series: Pandas Series of dispersions with gene identifiers.
    """
    counts = dds.X
    means = counts.mean(axis=0)
    variances = counts.var(axis=0)
    dispersions = variances / means
    dds.var['dispersion'] = dispersions
    return dds.var['dispersion']

def dispersion_trend(mean, a, b):
    """
    Model function for fitting the dispersion trend.
    
    Parameters:
    mean (float): Mean expression level.
    a (float): Coefficient a.
    b (float): Coefficient b.
    
    Returns:
    float: Fitted dispersion value.
    """
    return a / mean + b

def fit_dispersion_trend(dds):
    """
    Fit the dispersion trend across genes.
    
    Parameters:
    dds (anndata.AnnData): AnnData object containing the counts matrix and necessary metadata.

    Returns:
    np.ndarray: Array of fitted dispersions.
    """
    means = dds.X.mean(axis=0)
    dispersions = dds.var['dispersion'].values
    
    # Filter out zero means to avoid division by zero
    non_zero_means = means[means > 0]
    non_zero_dispersions = dispersions[means > 0]
    
    # Fit the model to the data
    popt, _ = curve_fit(dispersion_trend, non_zero_means, non_zero_dispersions)
    
    # Apply the fitted model to all means
    fitted_dispersions = dispersion_trend(means, *popt)
    dds.var['fitted_dispersion'] = fitted_dispersions
    
    return fitted_dispersions

def shrink_dispersions(dds):
    """
    Shrink the dispersions towards the fitted trend.
    
    Parameters:
    dds (anndata.AnnData): AnnData object containing the counts matrix and necessary metadata.

    Returns:
    pd.Series: Pandas Series of shrunken dispersions with gene identifiers.
    """
    raw_dispersions = dds.var['dispersion']
    fitted_dispersions = dds.var['fitted_dispersion']
    
    # Shrinkage formula: weighted average of raw and fitted dispersions
    shrunken_dispersions = 1 / (1 / raw_dispersions + 1 / fitted_dispersions)
    dds.var['shrunken_dispersion'] = shrunken_dispersions
    
    return dds.var['shrunken_dispersion']
