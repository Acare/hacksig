#' Single Sample Z-score
#'
#' Obtain single sample z-scores, as implemented in XXX, from a cohort of samples.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param signature A list of named character vectors corresponding to gene signatures.
#'
#' @return A tibble with the first column indicating sample identifiers (`sample_id`)
#'  and all the other columns ....
#' @export
#'
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
#' compute_zscore(test_expr, test_list)

compute_zscore <- function(expr_data, signatures) {

    single_sig_zscore <- function(dataset, genes) {
        filt_data <- dataset[genes, ]
        scaled_data <- base::scale(t(filt_data), center = TRUE, scale = TRUE)
        zscore_vec <- base::rowSums(scaled_data) / sqrt(length(genes))
        tibble::enframe(zscore_vec, name = "sample_id", value = "zscore")
    }

    result <- base::lapply(
        signatures,
        FUN = function(geneset) {
            single_sig_zscore(dataset = expr_data, genes = geneset)
        }
    )
    result <- dplyr::bind_rows(result, .id = "signature_id")
    tidyr::pivot_wider(result,
                       id_cols = sample_id,
                       names_from = signature_id,
                       values_from = zscore)

}
