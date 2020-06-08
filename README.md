# Deconvolution Algorithm Benchmarking

## (Basic) Input:  
* A list of wrapper around deconvolution models / algorithms
* scRNA-Seq data set (count matrix, and pheno information)
* RNA-Seq data with FACS

## Quality metrics: 
### On RNA-Seq data
* Correlation per cell type + overall correlation + bootstrapped correlation

### scRNA-Seq data  
* Missing reference profiles
* highly similar/nearly indistinguishable reference profiles
* wrong reference profiles in X, or the mixtures.  
* algorithms behaviour towards less genes, or less biological variance in the training data (For this, we have to retrain every algorithm.)  

## Output: 
RMD knitted html file. Or an interactive shiny visualization

# Aim of the package: 
## One model does not deconvolute every data best
I think that an deconvolution model that is fitted to a tissue scenario leads to the best results.  
However, the same model will underperform in another scenario. 
## DTD is the best algorithm to fit a decnvolution model
DTD sees artificial mixtures during training, and adapts the residuals accorindgly.  
Therefore, the benchmark package will ~~hopefully~~ show that DTD outperforms the other algorithm. 


# Include custom algorithms
## write a wrapper function for the new algorithm
The wrapper needs to support the following arguments:
  - `exprs`: non-negative numeric matrix containing single cell profiles as columns and features as rows.
  - `param`: pheno data.frame. Every row is a single cell pheno information and it can be assumed that `cell_type` is contained as column in the dataframe.
  - `bulks`: matrix containing bulk expression profiles as columns

`exprs` and `param` is single cell data that can be considered training data or from which the signature matrix should be inferred.
  
The function should return a list with two entries, `est.props` (the estimated proportions) and `sig.matrix`, the effectively used signature matrix.

## test your wrapper
call `DeconvolutionAlgorithmBenchmarking::check_algorithm(run_<ALGO>)`???

## use it in the benchmark
call the benchmark with a non-null `input.algorithms` parameter and include your wrapper as a list, such as `list(name="ALGO", algorithm=run_<ALGO>)` (where `run_<ALGO>` is your newly written wrapper function).

# Available benchmarks
All benchmarks get as parameters:
  - training data
  - test data
  - bulk data with its true proportions to score the performance of every algorithm
  - a set of algorithms to test
  - `n.repeats`, the number of repetitions for every algorithm (runs differ by their signature matrix and training data)
  - additional, benchmark specific parameters

## geneset benchmark
enabled, if `simulation.genes = TRUE` is passed to `benchmark`.

Compares the performance of the selected algorithms across different gene set definitions.

### detailed description
Expects a list of named genesets and runs every algorithm on the intersection of the genes available in the single cell data and the genes that are contained in the given gene sets.

## sample benchmark
enabled, if `simulation.samples = TRUE` is passed to `benchmark`.

Compares the influence of the size of the trainings set on the deconvolution performance across available algorithms.

### detailed description
Scans through fractions of the available single cell data to see how algorithms perform when the amount of training data (used for the determination of the signature matrix and potential training phases) is limited. The step width is given by `step.size` and the fraction `f` of selected cell profiles varies from `step.size` to its largest multiple smaller than 1. In every step, a fraction `f` of cells of each type are selected without replacement and added to the set of selected cells. Then, all algorithms are called with that set of cells and the deconvolution performance is determined on the bulks.

## bulk benchmark
enabled, if `simulation.bulks = TRUE` is passed to `benchmark`.

This "most basic" benchmark compares deconvolution performance of all selected algorithms by scoring how well the algorithms deconvolute given (real or simulated) bulks. See bulk simulation for details on how the bulks are simulated.

## subtype benchmark
enabled, if `simulation.subtypes = TRUE` is passed to `benchmark`.

The deconvolution performance is measured on simulated, "unbalanced" bulks. That means, cell populations of every individual type are not sampled uniformly but the sampled cell populations in the bulk are drawn from a different (sub-) population than were used for the signature matrix (full population). This may happen similarly in reality if certain subtypes dominate in the bulk but the coarse type should be determined, or likewise, if bulks were obtained from cells in a slightly different state than the single cell data (technical, batch or biological effects).

### detailed description
First, a joint tSNE embedding is computed on all single cell data and then performs k-means clustering on every cell type population in that (two-dimensional) space. By default, k = 3. Every cell profile is then assigned its cluster as a subtype. 


- TODO: bulk data does not "see" subtypes??? signature does?
- TODO: Simulated bulk has (on average!) same distribution of subtypes as sc input!


# Simulation of bulks
- TODO: I personally would find it better if not `fraction.per.bulk` was given but `ncells.per.bulk` (a typical bulk is a certain number of cells rather than a certain fraction of otherwise unrelated sc measurements)
- TODO: may include other ways of simulating bulks???

Every simulated bulk is averaged over a given fraction `fraction.per.bulk` of all single cells, i.e. the distribution of the single cell data determines the average distribution of the bulks. If `sum.to.count = TRUE`, every bulk profile is normed so that the sum over the gene expression is equal to the number of genes.
