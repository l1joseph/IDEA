from .dataset import IDEADataSet
from .differential_expression import IDEAResults, create_design_matrix, fit_negative_binomial_model, perform_wald_test, compute_log2_fold_change, runidea
from .dispersion import estimate_dispersions, fit_dispersion_trend, shrink_dispersions
from .normalization import median_of_ratios, size_factors, normalize_counts
from .plots import map_ensembl_to_gene_name, plot_volcano, plot_heatmap

__version__ = '0.1.0'
__all__ = [
    'IDEADataSet',
    'IDEAResults',
    'create_design_matrix', 
    'fit_negative_binomial_model',
    'perform_wald_test', 
    'compute_log2_fold_change', 
    'runidea',
    'estimate_dispersions',
    'fit_dispersion_trend',
    'shrink_dispersions',
    'median_of_ratios',
    'size_factors',
    'normalize_counts',
    'map_ensembl_to_gene_name',
    'plot_volcano',
    'plot_heatmap'
]
