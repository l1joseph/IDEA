# IDEA

## IDEA (Interactive Differential Expression Analysis) for CSE 185 Final Project

Gavin Simmons [gasimmons@ucsd.edu](mailto:gasimmons@ucsd.edu), Leo Joseph [l1joseph@ucsd.edu](mailto:l1joseph@ucsd.edu), Gayathri Donepudi [gdonepudi@ucsd.edu](mailto:gdonepudi@ucsd.edu)

## Project Description

IDEA (Interactive Differential Expression Analysis) is a Python package designed for performing differential expression analysis on gene expression data. The goal of IDEA is to provide researchers with a user-friendly and efficient tool to identify differentially expressed genes between different conditions or groups.

IDEA aims to be a Python equivalent to the widely used DESeq2 tool, offering similar functionality and performance. The package is implemented using Python, leveraging popular libraries such as NumPy, Pandas, and SciPy for data manipulation and statistical analysis.

## Installation

To install IDEA, you need to have Python 3.11 installed on your system. You can install IDEA using pip:

```
pip install git+https://github.com/l1joseph/IDEA.git
```

To run the example notebooks as well as development of this package, you will need a conda environment specific to IDEA. We have provided a yaml file with all the necessary dependencies installed. Run the below commands to activate it.

```
conda env create -f IDEA.yaml
conda activate IDEA
```

IDEA requires the following dependencies:

- AnnData
- Matplotlib
- NumPy
- Pandas
- SciPy
- Seaborn
- Statsmodels
- System

These dependencies will be automatically installed when you install IDEA using pip.

## Usage

To use IDEA, follow these steps:

1. Import the IDEA package in your Python script:

   ```python
   import idea
   ```

2. Load your gene expression data into an IDEADataSet object:

   ```python
   data = idea.IDEADataSet(counts, sample_info)
   ```

3. Convert data into an AnnData object:

   ```python
   adata = data.to_anndata()
   ```

4. Normalize data:

   ```python
   idea.size_factors(adata)
   idea.normalize_counts(adata)
   ```

5. Estimate dispersions:

   ```python
   idea.estimate_dispersions(adata)
   idea.fit_dispersion_trend(adata)
   idea.shrink_dispersions(adata)
   ```

6. Perform differential expression analysis:

   ```python
   contrast = ('condition', 'A', 'B')
   results = idea.runidea(adata, contrast)
   ```

7. Explore the results and visualize the differentially expressed genes:
   ```python
   idea_results = results.get_results()
   idea_results_mapped = idea.map_ensembl_to_gene_name(idea_results, 'path_to_file_with_gene_names')
   idea.plot_volcano(idea_results_mapped)
   idea.plot_heatmap(idea_results_mapped, adata)
   ```

## Examples

We provide a set of example Jupyter notebooks in the `example_notebook` prefix of the notebook name. These notebooks demonstrate how to use IDEA for various analysis tasks and provide step-by-step guides.
Current example notebooks:

- [example_notebook_Lab4Data.ipynb](example_notebook_Lab4Data.ipynb) (Data Obtained from [High fat diet-induced changes of mouse hepatic transcription and enhancer activity can be reversed by subsequent weight loss *Scientific Reports* 2017](https://www.nature.com/articles/srep40220.pdf))


## Contact

If you have any questions, suggestions, or issues, please feel free to contact us:

- Gavin Simmons [gasimmons](mailto:gasimmons@ucsd.edu)
- Leo Joseph [l1joseph@ucsd.edu](mailto:l1joseph@ucsd.edu)
- Gayathri Donepudi [gdonepudi@ucsd.edu](mailto:gdonepudi@ucsd.edu)

You can also open an issue on the project's GitHub repository.


## Future Plans

We have exciting plans for the future development of IDEA, including:

- Integration with popular visualization libraries for enhanced data exploration.
- Support for additional input data formats.
- Improved performance and scalability for large datasets.
- Integration with other bioinformatics tools and pipelines.

Stay tuned for upcoming releases and updates!

---
