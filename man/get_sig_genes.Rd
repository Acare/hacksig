% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_genes.R
\name{get_sig_genes}
\alias{get_sig_genes}
\title{Get signature gene identifiers}
\usage{
get_sig_genes(keywords = "all")
}
\arguments{
\item{keywords}{A character vector indicating keywords for a group of signatures.
The default (\code{"all"}) will cause the function to check for all the signatures
implemented in \code{hacksig}.}
}
\value{
A named list of gene signatures.
}
\description{
Obtain gene signatures implemented in \code{hacksig} as a named list of gene symbols.
}
\examples{
get_sig_genes()
get_sig_genes("estimate")
}
\seealso{
\code{\link[=get_sig_info]{get_sig_info()}} to get valid keywords for signatures.
}
