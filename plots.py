import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plotVolcano(res):
    # Create a volcano plot
    plt.figure(figsize=(10, 8))
    res['-log10(p_value)'] = -np.log10(res['p_value'])
    significant = (res['p_value'] < 0.05) & (abs(res['log2_fold_change']) > 1)
    
    plt.scatter(res[~significant]['log2_fold_change'], res[~significant]['-log10(p_value)'], color='gray', alpha=0.5)
    plt.scatter(res[significant]['log2_fold_change'], res[significant]['-log10(p_value)'], color='red', alpha=0.8)
    
    plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--', linewidth=0.8)
    plt.axvline(x=1, color='blue', linestyle='--', linewidth=0.8)
    plt.axvline(x=-1, color='blue', linestyle='--', linewidth=0.8)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(p-value)')
    plt.title('Volcano Plot')
    plt.show()

def plotHeatmap(res, adata, top_n=10):
    # Select top N differentially expressed genes based on p-value
    top_genes = res.nsmallest(top_n, 'p_value')['gene']
    counts = pd.DataFrame(adata.layers['normalized'], index=adata.obs.index, columns=adata.var.index)
    top_counts = counts[top_genes]

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(top_counts.T, cmap='viridis', xticklabels=adata.obs['condition'], yticklabels=top_genes)
    plt.title(f'Top {top_n} Differentially Expressed Genes')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.show()
