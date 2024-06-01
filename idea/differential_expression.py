from typing import Tuple

import anndata
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.stats.multitest as smm


class IDEAResults:
    def __init__(self, dds, contrast):
        self.dds = dds
        self.contrast = contrast
        self.results = None
    
    def get_results(self):
        return self.results
        

def create_design_matrix(phenotype_data: pd.DataFrame, condition: str) -> pd.DataFrame:
    """
    Create design matrix with dummy variables for the specified condition.

    Parameters
    ----------
    phenotype_data : pd.DataFrame
        DataFrame containing the phenotype data.
    condition : str
        Name of the column containing the condition information.

    Returns
    -------
    pd.DataFrame
        Design matrix with dummy variables for the specified condition.
    """
    return pd.get_dummies(phenotype_data[condition], drop_first=False).astype(np.int64)

def fit_negative_binomial_model(counts: np.ndarray, design_matrix: pd.DataFrame) -> sm.GLM:
    """
    Fit a negative binomial regression model for the given counts and design matrix.

    Parameters
    ----------
    counts : np.ndarray
        Array of count data.
    design_matrix : pd.DataFrame
        Design matrix for the regression model.

    Returns
    -------
    sm.GLM
        Fitted negative binomial regression model.
    """
    return sm.GLM(counts, design_matrix, family=sm.families.NegativeBinomial()).fit()

def perform_wald_test(model: sm.GLM, contrast_matrix: np.ndarray) -> Tuple[float, float, float]:
    """
    Perform a Wald test using the specified contrast matrix.

    Parameters
    ----------
    model : sm.GLM
        Fitted negative binomial regression model.
    contrast_matrix : np.ndarray
        Contrast matrix for the Wald test.

    Returns
    -------
    Tuple[float, float, float]
        Tuple containing the Wald test statistic, p-value, and p-value (unadjusted).
    """
    wald_test = model.wald_test(contrast_matrix)
    if wald_test is None or not hasattr(wald_test, 'statistic'):
        return np.nan, np.nan, np.nan

    statistic = wald_test.statistic.item()
    p_value = wald_test.pvalue.item()
    return statistic, p_value, p_value

def compute_log2_fold_change(model: sm.GLM, design_matrix: pd.DataFrame, level1: str, level2: str) -> float:
    """
    Compute the log2 fold change between two levels of a condition.

    Parameters
    ----------
    model : sm.GLM
        Fitted negative binomial regression model.
    design_matrix : pd.DataFrame
        Design matrix for the regression model.
    level1 : str
        Name of the first level of the condition.
    level2 : str
        Name of the second level of the condition.

    Returns
    -------
    float
        Log2 fold change between the two levels of the condition.
    """
    level1_coef = model.params[design_matrix.columns.get_loc(level1)]
    level2_coef = model.params[design_matrix.columns.get_loc(level2)]
    return level2_coef - level1_coef

def runidea(adata: anndata.AnnData, contrast: Tuple[str, str, str]) -> pd.DataFrame:
    """
    Perform differential expression analysis using negative binomial regression.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing the counts matrix and necessary metadata.
    contrast : Tuple[str, str, str]
        Tuple specifying the contrast for differential expression (e.g., ('condition', 'A', 'B')).

    Returns
    -------
    pd.DataFrame
        DataFrame containing the results of the differential expression analysis.
    """
    counts = adata.layers['normalized']
    phenotype_data = adata.obs
    genes = adata.var.index
    condition, level1, level2 = contrast

    design_matrix = create_design_matrix(phenotype_data, condition)

    results = []
    p_values = []
    for i, gene in enumerate(genes):
        y = counts[:, i]
        model = fit_negative_binomial_model(y, design_matrix)

        contrast_matrix = np.zeros(len(design_matrix.columns), dtype=np.int64)
        contrast_matrix[design_matrix.columns.get_loc(level1)] = -1
        contrast_matrix[design_matrix.columns.get_loc(level2)] = 1

        statistic, p_value, _ = perform_wald_test(model, contrast_matrix)
        log2_fold_change = compute_log2_fold_change(model, design_matrix, level1, level2)

        if not np.isnan(statistic) and not np.isnan(p_value):
            results.append({'gene': gene, 'log2_fold_change': log2_fold_change,
                            'statistic': statistic, 'p_value': p_value})
            p_values.append(p_value)
        else:
            logger.warning(f"No valid Wald test result for gene: {gene}")

    # Remove NaN values from p_values before multiple testing correction
    p_values = [p for p in p_values if not np.isnan(p)]
    _, adjusted_p_values, _, _ = smm.multipletests(p_values, method='fdr_bh')

    results_df = pd.DataFrame(results)
    results_df['padj'] = adjusted_p_values

    # Create IDEAResults object
    idea_results = IDEAResults(adata, contrast)
    idea_results.results = results_df

    return idea_results
