import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def map_ensembl_to_gene_name(de_results, gene_names_file):
    """
    Map Ensembl IDs to gene names and remove genes with no mapping.

    Parameters
    ----------
    de_results : pd.DataFrame
        DataFrame containing the differential expression analysis results.
    gene_names_file : str
        Path to the file containing the mapping between Ensembl IDs and gene names.

    Returns
    -------
    pd.DataFrame
        DataFrame with Ensembl IDs mapped to gene names and genes with no mapping removed.
    """
    # Read the gene names file
    gene_names = pd.read_csv(gene_names_file, sep='\t', header=None, names=['ensembl_id', 'gene_name'])

    # Create a dictionary for mapping Ensembl IDs to gene names
    id_to_name = dict(zip(gene_names['ensembl_id'], gene_names['gene_name']))

    # Map Ensembl IDs to gene names in the results DataFrame
    de_results['gene_name'] = de_results['gene'].map(id_to_name)

    # Remove genes with no mapping
    de_results = de_results.dropna(subset=['gene_name'])

    return de_results

def plot_volcano(results_df, log2_fc_col='log2_fold_change', pval_col='p_value', alpha=0.05, lfc_threshold=1.0):
    """
    Generate a volcano plot for the differential expression results.
    Parameters:
    results_df (pd.DataFrame): DataFrame containing the differential expression analysis results.
    log2_fc_col (str): Column name for log2 fold change.
    pval_col (str): Column name for p-value.
    alpha (float): Significance threshold for p-value.
    lfc_threshold (float): Threshold for log2 fold change to consider a gene significantly differentially expressed.
    """
    # Ensure p-value column is numeric
    results_df[pval_col] = pd.to_numeric(results_df[pval_col], errors='coerce')
    # Calculate -log10(p-value)
    results_df['-log10_pvalue'] = -np.log10(results_df[pval_col])
    # Determine significance and direction of regulation
    results_df['significance'] = (results_df[pval_col] < alpha) & (np.abs(results_df[log2_fc_col]) > lfc_threshold)
    results_df['regulation'] = ['upregulated' if lfc > 0 and sig else 'downregulated' if lfc < 0 and sig else 'nonsignificant'
                                for lfc, sig in zip(results_df[log2_fc_col], results_df['significance'])]
    
    # Get the top significantly expressed genes
    top_sig_genes = results_df[results_df['significance']].nsmallest(10, 'padj')['gene_name'].tolist()
    if len(top_sig_genes) < 10:
        top_sig_gene_labels = top_sig_genes
    else:
        top_sig_gene_labels = [f'{gene_name[:15]}...' for gene_name in top_sig_genes]
    
    # Create the volcano plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=log2_fc_col, y='-log10_pvalue', data=results_df,
                    hue='regulation', palette={'upregulated': 'red', 'downregulated': 'blue', 'nonsignificant': 'gray'},
                    alpha=0.6)
    plt.axhline(-np.log10(alpha), ls='--', color='black', lw=0.5)
    plt.axvline(lfc_threshold, ls='--', color='black', lw=0.5)
    plt.axvline(-lfc_threshold, ls='--', color='black', lw=0.5)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 p-value')
    plt.title('Volcano Plot of Differential Expression')
    plt.legend(title='Regulation', loc='upper right')
    
    # Add gene names for top significantly expressed genes
    for i, gene_name in enumerate(top_sig_gene_labels):
        plt.annotate(gene_name, (results_df.loc[results_df['gene_name'] == top_sig_genes[i], log2_fc_col].values[0],
                                 results_df.loc[results_df['gene_name'] == top_sig_genes[i], '-log10_pvalue'].values[0]),
                     fontsize=8)
    
    plt.show()

def plot_heatmap(results_df, adata, top_n=20):
    """
    Plot a heatmap of the top differentially expressed genes.
    Parameters:
    results_df (pd.DataFrame): DataFrame containing the differential expression analysis results.
    adata (anndata.AnnData): AnnData object containing the counts matrix and metadata.
    top_n (int): Number of top genes to display in the heatmap.
    """
    # Select the top_n differentially expressed genes by adjusted p-value
    top_genes = results_df.nsmallest(top_n, 'padj')['gene']
    top_genes_data = adata[:, top_genes].X
    
    # Ensure the data is in the correct format
    if not isinstance(top_genes_data, np.ndarray):
        top_genes_data = top_genes_data.toarray()
    
    # Get the corresponding gene names
    top_gene_names = results_df.loc[results_df['gene'].isin(top_genes), 'gene_name'].tolist()
    
    # Create a heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(top_genes_data, yticklabels=adata.obs.index, xticklabels=top_gene_names, cmap='RdBu_r', cbar=True)
    plt.title(f'Top {top_n} Differentially Expressed Genes')
    plt.xlabel('Gene Names')
    plt.ylabel('Samples')
    plt.show()
