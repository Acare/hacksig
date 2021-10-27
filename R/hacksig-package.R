#' hacksig: A Tidy Framework to Hack Gene Expression Signatures in R
#'
#' @description
#' The hacksig package has been designed for the purpose of simplifying the way
#' in which gene expression signatures scores are computed. It makes use of three
#' different single sample score calculation methods, namely z-score, single sample
#' GSEA and singscore. It is a manually curated collection of gene expression
#' signatures found in literature, such as CINSARC, ESTIMATE and Immunophenoscore.
#' Moreover, it supports parallelization through the `future` framework.
#'
#' @section Get gene signatures scores in different ways:
#' The main function of the package is \code{\link{hack_sig}} and it can be used to:
#'   * obtain single sample scores with one of three methods (z-score, ssGSEA,
#'     singscore) for a custom list of gene signatures;
#'   * obtain single sample scores for literature gene signatures implemented in
#'     the package, either with the default method or with one of the three single
#'     sample methods.
#'
#' @section Check if gene signatures are applicable:
#' Sometimes your gene expression matrix can miss some genes due to some prior
#' filtering procedure. The function \code{\link{check_sig}} can be used to check
#' the number of genes present in your gene expression matrix either for your custom
#' gene signatures or for those implemented in `hacksig`.
#'
#' @author Andrea Carenzo, Loris De Cecco, Federico Pistore
#' @docType package
#' @name hacksig
#' @aliases hacksig-package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
