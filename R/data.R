#' Gene-level information for each signature implemented
#'
#' A data set used by the `hacksig` functions to retrieve information about
#' signatures implemented in the package.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{signature_id}{a unique ID given to the signature}
#'   \item{gene_symbol}{Gene Symbol ID}
#'   \item{entrez_gene_id}{Entrez gene ID}
#'   \item{gene_type}{gives the gene group}
#'   \item{gene_weight}{gives the gene weight in the model}
#'   \item{method}{default method implemented in the original work}
#' }
"signatures_data"


#' A toy gene expression matrix
#'
#' A gene expression matrix ...
#'
#' @format A matrix with 20000 random gene symbols as rows and 20 samples as columns.
"test_expr"

#' Hallmark gene set list
#'
#' A gene expression matrix ...
#'
#' @format A named list of 50 signatures corresponding to the Hallmark gene sets.
"hallmark_list"
