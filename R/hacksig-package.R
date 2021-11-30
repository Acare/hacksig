#' hacksig: A Tidy Framework to Hack Gene Expression Signatures
#'
#' @description
#' The `hacksig` package has been designed for the purpose of simplifying the way
#'   in which gene expression signatures scores are computed. It is a manually
#'   curated collection of gene expression signatures found in literature and makes
#'   use of three different single sample score calculation methods. Moreover,
#'   parallel computation is supported through the `future` framework.
#' @section Get gene signatures scores in different ways:
#' The main function of the package is [hack_sig()] and it can be used to:
#'
#'   * obtain single sample scores with one of three methods (z-score, ssGSEA,
#'     singscore) for a custom list of gene signatures;
#'   * obtain single sample scores for a number of manually curated gene signatures
#'     either with the original publication method or with one of the three single
#'     sample methods.
#'
#'   Once single sample scores are obtained, you can assign your samples into
#'   signature classes with [hack_class()].
#'
#'   In addition, other more complex methods are implemented through:
#'
#'   * [hack_cinsarc()], for the CINSARC classification;
#'   * [hack_estimate()], for the ESTIMATE method;
#'   * [hack_immunophenoscore()], for the Immunophenoscore.
#'
#'   Information about implemented signatures can be obtained with [get_sig_info()].
#' @section Check if gene signatures are applicable to your data:
#'   Sometimes your gene expression matrix can miss some genes due to some prior
#'   filtering procedure. The function [check_sig()] can be used to check how many
#'   genes your expression matrix contain for every input signature.
#' @author Andrea Carenzo, Loris De Cecco, Federico Pistore
#' @docType package
#' @name hacksig
#' @aliases hacksig-package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
