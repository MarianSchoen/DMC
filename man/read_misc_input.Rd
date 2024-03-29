% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h5_read_misc.R
\name{read_misc_input}
\alias{read_misc_input}
\title{read input to benchmark function from hdf5 file}
\usage{
read_misc_input(filename)
}
\arguments{
\item{filename}{string, which file should be read in?}
}
\value{
list with:
   - 'genesets': list of string vectors
   - 'algorithms': vector containing algorithm names
   - 'grouping': numeric vector grouping samples in test and training set
   - 'function.call': call object
}
\description{
read input to benchmark function from hdf5 file
}
