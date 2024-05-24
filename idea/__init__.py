from .dataset import IDEADataSet
from .differential_expression import IDEAResults, deseq, adjust_pvalues
from .dispersion import estimate_dispersions, fit_dispersion_trend, shrink_dispersions
from .normalization import median_of_ratios, size_factors, normalize_counts
from .plots import plotVolcano, plotHeatmap

__version__ = '0.1.0'
__all__ = [
    'IDEADataSet',
    'IDEAResults',
    'deseq',
    'adjust_pvalues',
    'estimate_dispersions',
    'fit_dispersion_trend',
    'shrink_dispersions',
    'median_of_ratios',
    'size_factors',
    'normalize_counts',
    'plotVolcano',
    'plotHeatmap'
]