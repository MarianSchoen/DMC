% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h5_read_data.R
\name{read_data}
\alias{read_data}
\title{read data conforming to the format used by \link{benchmark} for data storage}
\usage{
read_data(filename)
}
\arguments{
\item{filename}{string, which file should be read in?}
}
\value{
list with: \cr
 - sc.counts: numeric matrix, features as rows,
   scRNA-Seq profiles as colums\cr
 - sc.pheno: dataframe with scRNA-Seq profiles as rows, and multiple pheno
  columns \cr
 - bulk.counts: numeric matrix, features as rows, bulk measurements as
 columns\cr
 - bulk.props: numeric matrix, cell types as rows, bulks as columns\cr
 - bulk.pheno: dataframe containing pheno data for bulks in columns,
 bulks as rows
}
\description{
read data conforming to the format used by \link{benchmark} for data storage
}
