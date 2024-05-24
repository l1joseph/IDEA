# IDEA

## IDEA (Interactive Differential Expression Analysis) for CSE 185 Final Project

Gavin Simmons [gasimmons@ucsd.edu](mailto:gasimmons@ucsd.edu), Leo Joseph [l1joseph@ucsd.edu](mailto:l1joseph@ucsd.edu), Gayathri Donepudi [gdonepudi@ucsd.edu](mailto:gdonepudi@ucsd.edu)

## Project Description

IDEA (Interactive Differential Expression Analysis) is a Python package designed for performing differential expression analysis on gene expression data. The goal of IDEA is to provide researchers with a user-friendly and efficient tool to identify differentially expressed genes between different conditions or groups.

IDEA aims to be a Python equivalent to the widely used DESeq2 tool, offering similar functionality and performance. The package is implemented using Python, leveraging popular libraries such as NumPy, Pandas, and SciPy for data manipulation and statistical analysis.

## Installation

To install IDEA, you need to have Python 3.11 installed on your system. You can install IDEA using pip:

```
#implement this after we finalize package details
# pip install idea

pip install git+https://github.com/l1joseph/IDEA.git
# or whatever we go with
```

To run the example notebooks as well as development of this package, you will need a conda environment specific to IDEA. We have prodived a yaml file with all the necessary dependencies installed. Run the below commands to activate it.

```
conda env create -f IDEA.yaml
conda activate IDEA
```

IDEA requires the following dependencies:

- NumPy
- Pandas
- SciPy
- (CONTINUE AS WE ADD MORE REQS)

These dependencies will be automatically installed when you install IDEA using pip.

## Usage

To use IDEA, follow these steps:

1. Import the IDEA package in your Python script:

   ```python
   import idea
   ```

2. Load your gene expression data into a Pandas DataFrame:

   ```python
   data = idea.load_data('expression_data.tsv')
   ```

3. Perform differential expression analysis:

   ```python
   results = idea.differential_expression(data, group1, group2)
   ```

4. Explore the results and visualize the differentially expressed genes:
   ```python
   idea.plot_volcano(results)
   idea.plot_heatmap(results)
   ```

For more detailed usage examples and documentation, please refer to the [IDEA Documentation](https://idea.readthedocs.io/).
(NO CLUE IF WE'RE GOING TO DO THIS)

## Examples

We provide a set of example Jupyter notebooks in the `examples/` directory of the project repository. These notebooks demonstrate how to use IDEA for various analysis tasks and provide step-by-step guides.

## Documentation

The complete documentation for IDEA can be found at [https://idea.readthedocs.io/](https://idea.readthedocs.io/). It includes API references, user guides, and tutorials.
(NO CLUE IF WE'RE GOING TO DO THIS)

## Contributing

We welcome contributions to IDEA! If you would like to contribute, please follow these guidelines:

1. Fork the repository and create a new branch for your feature or bug fix.
2. Write clear and concise commit messages.
3. Submit a pull request describing your changes and their purpose.

For more detailed information, please refer to our [Contributing Guide](CONTRIBUTING.md).

## License

(I HAVEN'T SET A LICENSE FOR THIS WHEN I INIT THE REPO, DON'T KNOW IF WE HAVE TO)
IDEA is distributed under the MIT License. See [LICENSE](LICENSE) for more information.

## Contact

If you have any questions, suggestions, or issues, please feel free to contact us:

- Gavin Simmons [gasimmons](mailto:gasimmons@ucsd.edu)
- Leo Joseph [l1joseph@ucsd.edu](mailto:l1joseph@ucsd.edu)
- Gayathri Donepudi [gdonepudi@ucsd.edu](mailto:gdonepudi@ucsd.edu)

You can also open an issue on the project's GitHub repository.

## Acknowledgements

(ADD ACKNOWLEDGEMENTS HERE)

## Future Plans

We have exciting plans for the future development of IDEA, including:

- Integration with popular visualization libraries for enhanced data exploration.
- Support for additional input data formats.
- Improved performance and scalability for large datasets.
- Integration with other bioinformatics tools and pipelines.

Stay tuned for upcoming releases and updates!

---
