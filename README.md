# Deconvolution Algorithm Benchmarking

## (Basic) Input:  
* A list of wrapper around deconovlution models
* scRNA-Seq data set (count matrix, and pheno information)
* RNA-Seq data with FACS

## Quality metrics: 
### On RNA-Seq data
* Correlation per cell type + overall correlation 

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

