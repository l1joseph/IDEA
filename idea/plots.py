import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plotVolcano(results_df, log2_fc_col='log2_fold_change', pval_col='p_value', alpha=0.05, lfc_threshold=1):
    """
    Generate a volcano plot for the differential expression results.
    
    Parameters:
    - results_df: DataFrame with differential expression results.
    - log2_fc_col: Column name for log2 fold change values.
    - pval_col: Column name for p-values.
    - alpha: Significance level for p-values.
    - lfc_threshold: Threshold for log2 fold change to consider a gene significantly differentially expressed.
    """
    results_df['-log10_pvalue'] = -np.log10(results_df[pval_col])
    
    # Determine upregulated and downregulated genes
    results_df['regulation'] = 'not_significant'
    results_df.loc[(results_df[log2_fc_col] > lfc_threshold) & (results_df[pval_col] < alpha), 'regulation'] = 'upregulated'
    results_df.loc[(results_df[log2_fc_col] < -lfc_threshold) & (results_df[pval_col] < alpha), 'regulation'] = 'downregulated'
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=results_df, x=log2_fc_col, y='-log10_pvalue', hue='regulation', palette={
        'upregulated': 'red', 'downregulated': 'blue', 'not_significant': 'grey'
    })
    
    plt.axhline(y=-np.log10(alpha), linestyle='--', color='black', linewidth=1)
    plt.axvline(x=lfc_threshold, linestyle='--', color='black', linewidth=1)
    plt.axvline(x=-lfc_threshold, linestyle='--', color='black', linewidth=1)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 p-value')
    plt.title('Volcano Plot')
    plt.legend(title='Regulation')
    plt.show()

def plotHeatmap(results_df, adata, num_genes=20, log2_fc_col='log2_fold_change'):
    """
    Generate a heatmap of the top differentially expressed genes.
    
    Parameters:
    - results_df: DataFrame with differential expression results.
    - adata: AnnData object with the normalized counts data.
    - num_genes: Number of top differentially expressed genes to include in the heatmap.
    - log2_fc_col: Column name for log2 fold change values.
    """
    # Select top differentially expressed genes based on log2 fold change
    top_genes = results_df.sort_values(by=log2_fc_col, ascending=False).head(num_genes).index
    selected_data = adata[:, top_genes].X
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(selected_data, cmap='vlag', xticklabels=top_genes, yticklabels=adata.obs.index)
    plt.xlabel('Genes')
    plt.ylabel('Samples')
    plt.title(f'Top {num_genes} Differentially Expressed Genes')
    plt.show()
