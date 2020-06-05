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
