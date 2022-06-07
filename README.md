# Deconvolution Model Comparison  

One model does not deconvolute all datasets best. There are several factors, biological and technical, that influence deconvolution performance, which may have different effects on different algorithms.  
This package enables easy comparison of deconvolution models on a given dataset, which may be used to determine the best algorithm for a specific use case.

## Contributors
* Marian Schön
* Tim Mirus
* Jakob Simeth

# Details

## Input:  
* A list of wrapper around deconvolution models / algorithms
* scRNA-Seq data set (count matrix, and pheno information)
* bulk RNA-Seq data with FACS
* ...

## Output: 
RMD knitted html file. Plots are saved to specified working directory under 'report_plots'

## Quality metrics: 
Deconvolution performance is determined as the pearson correlation coefficient between the real and estimated cell type quantities for a given cell type across all bulks in the dataset:  
```math
r_{ct} = cor(C_{ct,.}, \hat{C}_{ct,.})
```
The total performance of an algorithm is defined as the mean performance across all celltypes.

## scRNA-Seq data 
When no RNA-Seq data set with ground-truth cell type quantities is available, a number of benchmarks can be performed based on single-cell data only:

* Performance on simulated bulks
* Performance on simulated bulks with varying amounts of training data (single-cell profiles)
* Performance on simulated bulks with different gene sets
* Performance on simulated bulks with finer cell type labels obtained by clustering of single-cell profiles

## Simulation of bulks

Every simulated bulk is averaged over a given fraction `fraction.per.bulk` of all single cells, i.e. the distribution of the single cell data determines the average distribution of the bulks. If `sum.to.count = TRUE`, every bulk profile is normed so that the sum over the gene expression is equal to the number of genes. The amount of cells of each type included in a certain bulks is distributed non-uniformly, because the created bulks would otherwise all reflect the cell type proportions of the single-cell data, which is not a realistic scenario.
    
# Usage

## Include custom algorithms

### Write a wrapper function for the new algorithm  

The wrapper needs to support the following arguments:
  - `exprs`: non-negative numeric matrix containing single cell profiles as columns and features as rows.
  - `param`: pheno data.frame. Every row is a single cell pheno information and it can be assumed that `cell_type` is contained as column in the dataframe.
  - `bulks`: matrix containing bulk expression profiles as columns

`exprs` and `param` is single cell data that can be considered training data or from which the signature matrix should be inferred.
The function must return a list with two entries, `est.props` (the estimated proportions) and `sig.matrix`, the effectively used signature matrix.

## Test your wrapper
The package contains a function for checking whether a wrapper function complies with the package standards
```
DMC::check_algorithm(list(algorithm = run_<ALGO>, name = "ALGO", model = NULL))
```

## Use it in the benchmark
Call the benchmark with the `input.algorithms` parameter and include your wrapper as a list, such as `list(name="ALGO", algorithm=run_<ALGO>)` (where `run_<ALGO>` is your newly written wrapper function).

# Available benchmarks
All benchmarks get as parameters:
  - training data
  - test data
  - bulk data with its true proportions to score the performance of every algorithm
  - a set of algorithms to test
  - the number of repetitions for every algorithm (runs differ by their signature matrix and training data)
  - additional, benchmark specific parameters

## Real data
A simple deconvolution of given bulk data, using all available single-cell data for training. The resulting cell type quantities are then compared to the ground-truth by correlation (see above).

## Bulk benchmark
* enabled if `simulation.bulks = TRUE` is passed to `benchmark`.

This is basically the same as the benchmark on real data, but with artificial bulks created from single-cell profiles. Useful if no RNA-Seq data with ground-truth cell type quantities is available.

## Geneset benchmark
* enabled if `simulation.genes = TRUE` is passed to `benchmark`.

Compares the performance of the selected algorithms across different gene set definitions / signatures.

## Sample benchmark
* enabled if `simulation.samples = TRUE` is passed to `benchmark`.

Compares the influence of the size and composition of the training set on the deconvolution performance across available algorithms by repeated random sub-sampling of the training set to different sizes.

## Subtype benchmark
* enabled if `simulation.subtypes = TRUE` is passed to `benchmark`.

Based on the given cell type labels, divide the single-cell profiles further into subtypes using hierarchical clustering based on a t-SNE embedding. This happens at different depths, creating cell type labels of increasing granularity. For each level of cell types, artificial bulks are deconvoluted to determine how well each algorithm can distinguish cell types of increasing similarity.

Further details on usage may be taken from the Vignette.
