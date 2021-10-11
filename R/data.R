#' Gene-level information for each signature implemented.
#'
#' A dataset containing information about signatures implemented in the `hacksig`
#' package.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{signature_id}{a unique ID given to the signature}
#'   \item{ensembl_gene_id}{ENSEMBL gene ID}
#'   \item{gene_symbol}{Gene Symbol ID}
#'   \item{entrez_gene_id}{Entrez gene ID}
#'   \item{gene_type}{gives the gene group}
#'   \item{gene_weight}{gives the gene weight in the model}
#'   \item{method}{default method implemented in the original work}
#' }
"signatures_data"
