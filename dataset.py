import pandas as pd
import numpy as np
import anndata as ad

class DESeqDataSet:
    def __init__(self, counts, sample_info):
        self.counts = counts
        self.sample_info = sample_info
        self.normalizedCounts = None
        self.sizeFactors = None
        self.meanVarZero = None
        self.designMatrix = None

    def get_counts(self):
        return self.counts

    def get_sample_info(self):
        return self.sample_info
    
    def to_anndata(self):
        # Create an AnnData object
        adata = ad.AnnData(X=self.counts.values, obs=self.sample_info, var=pd.DataFrame(index=self.counts.columns))
        return adata
