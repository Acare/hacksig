#' Singscore
#'
#' Obtain single sample GSEA scores for a list of gene signatures,
#' as implemented in *Foroutan et al., 2018*.
#'
#' @inheritParams compute_zscore
#' @return A tibble with the first column indicating sample identifiers (`sample_id`)
#'  and all the other columns ....
#' @references
#' Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., &
#' Davis, M. J. (2018). Single sample scoring of molecular phenotypes.
#' *BMC bioinformatics*, 19(1), 404.
#' [doi: 10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4).
#'
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation
#' analysis for microarray and RNA-seq data. *BMC bioinformatics*, 14, 7.
#' [doi: 10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7).
#'
#' @seealso \code{\link[singscore]{simpleScore}}
#' @examples
#' set.seed(123)
#'
#' test_expr <- data.frame(
#'     sample1 = rnorm(10, mean = 8, sd = 2),
#'     sample2 = rnorm(10, mean = 7, sd = 2),
#'     sample3 = rnorm(10, mean = 8, sd = 3),
#'     row.names = paste0("gene", 1:10)
#' )
#'
#' test_list <- list(
#'     gs1 = paste0("gene", c(1, 2, 6)),
#'     gs2 = paste0("gene", 1:5),
#'     gs3 = paste0("gene", c(3, 1, 9)),
#'     gs4 = paste0("gene", c(1, 3, 9))
#' )
#'
#' compute_singscore(test_expr, test_list)
#'
#' @export
compute_singscore <- function(expr_data, signatures) {
    single_sig_singscore <- function(dataset, genes) {
        rank_data <- apply(dataset, MARGIN = 2,
                           FUN = rank, na.last = "keep", ties.method = "average")
        singscore_vec <- apply(rank_data, MARGIN = 2,
                            FUN = es_ssgsea, geneset = genes)
        tibble::enframe(singscore_vec, name = "sample_id", value = "singscore")
    }
    result <- base::lapply(
        signatures,
        FUN = function(genes) {
            single_sig_singscore(dataset = expr_data, genes = genes)
        }
    )
    tidyr::pivot_wider(result,
                       id_cols = sample_id,
                       names_from = signature_id,
                       values_from = singscore)
}
