# scPositioner v1.1.0

## Global optimal-based spatial mapping of single-cell multi-omics

[![python >=3.9](https://img.shields.io/badge/python-%3E%3D3.9-brightgreen)](https://www.python.org/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13161096.svg)](https://doi.org/10.5281/zenodo.13161096)

scPositioner is a computational method to map single cells into spatial context and integrate multi omics

![avatar](images/workflow.jpg)


### Create and activate conda environment with requirements installed.
For scPositioner, the Python version need is over 3.9. If you have already installed a lower version of Python, consider installing Anaconda, and then you can create a new environment.
```
cd scPositioner-main

conda env create -f environment.yml -n scpositioner
conda activate scpositioner
```


### Install scPositioner
```
python setup.py build
python setup.py install
```

## Tutorials (single cell spatial mapping)

scPositioner requires the single-cell spatial omics data (stored as `.h5ad` format) as input, where cell population label of each cell needs to be provided. For spot-level SRT datasets, the deconvolution step should be done in advance.

Here is an example of scPositioner on spot-level SRT reference (10X Visium):
* [Demonstration of scPositioner on the spot-level data](tutorial/tutorial_kidney.ipynb)

An example of scPositioner on single-cell SRT reference (10X Xenium):
* [Demonstration of scPositioner on the single-cell data](tutorial/tutorial_BC.ipynb)

And an example of scPositioner on multi-omics (snRNA-seq, snATAC-seq, 10X Visium):
* [Demonstration of scPositioner on the multi-omics data](tutorial/tutorial_heart.ipynb)


## Acknowledgements


## About
scPositioner is developed by Jingyang Qian and Hudong Bao. Should you have any questions, please contact Hudong Bao at baohd@zju.edu.cn.

