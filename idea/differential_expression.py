import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import rankdata

class IDEAResults:
    def __init__(self, dds, contrast):
        self.dds = dds
        self.contrast = contrast
        self.results = None
    
    def get_results(self):
        return self.results
        
def deseq(dds, contrast):
    """
    Perform differential expression analysis.

    Parameters:
    dds (anndata.AnnData): AnnData object containing the counts matrix and necessary metadata.
    contrast (tuple): A tuple specifying the contrast for differential expression (e.g., ('condition', 'A', 'B')).

    Returns:
    DESeqResults: Object containing the results of the differential expression analysis.
    """
    # Extract counts and metadata
    counts = dds.layers['normalized']
    phenotype_data = dds.obs
    genes = dds.var.index
    condition = contrast[0]
    level1 = contrast[1]
    level2 = contrast[2]

    # Debugging output to ensure correct levels in phenotype data
    print("Phenotype data:")
    print(phenotype_data[condition])
    print(f"Expected levels: {level1}, {level2}")

    # Create design matrix
    design_matrix = pd.get_dummies(phenotype_data[condition], drop_first=False)
    print("Design matrix:")
    print(design_matrix.head())

    # Check if contrast levels are present in design matrix
    if level1 not in design_matrix.columns or level2 not in design_matrix.columns:
        raise ValueError(f"Contrast levels {level1} and {level2} not found in phenotype data.")

    # Recode contrast levels to 0 and 1 for analysis
    design_matrix = design_matrix[[level1, level2]].astype(int)
    design_matrix.columns = [level1, level2]
    print("Recoded design matrix:")
    print(design_matrix.head())

    # Fit the model for each gene
    results = []
    for i, gene in enumerate(genes):
        y = counts[:, i]
        model = sm.OLS(y, design_matrix).fit()
        summary = model.summary2().tables[1]
        log2_fold_change = summary.loc[level2, 'Coef.']
        p_value = summary.loc[level2, 'P>|t|']
        results.append({
            'gene': gene,
            'log2_fold_change': log2_fold_change,
            'p_value': p_value
        })

    results_df = pd.DataFrame(results)
    results_df['padj'] = adjust_pvalues(results_df['p_value'])

    # Create DESeqResults object
    deseq_results = DESeqResults(dds, contrast)
    deseq_results.results = results_df

    return deseq_results

def adjust_pvalues(pvalues):
    """
    Adjust p-values using the Benjamini-Hochberg method.

    Parameters:
    pvalues (pd.Series): Series of p-values.

    Returns:
    pd.Series: Series of adjusted p-values.
    """
    pvalues_sorted = np.sort(pvalues)
    pvalues_order = np.argsort(pvalues)
    n = len(pvalues)
    adjusted_pvalues = np.zeros(n)
    cumulative_min = np.inf
    for i in range(n-1, -1, -1):
        pval = pvalues_sorted[i]
        cumulative_min = min(cumulative_min, pval * n / (i + 1))
        adjusted_pvalues[pvalues_order[i]] = cumulative_min
    return adjusted_pvalues
