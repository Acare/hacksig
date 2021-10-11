#' Single sample z-score
#'
#' Obtain single sample z-scores for a list of gene signatures,
#' as implemented in *Lee et al. (2008)*.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param signatures A list of named character vectors corresponding to gene signatures.
#'
#' @return A tibble with the first column indicating sample identifiers (`sample_id`)
#'  and all the other columns ....
#'
#' @references
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation
#' analysis for microarray and RNA-seq data. *BMC bioinformatics*, 14, 7.
#' [doi :10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7).
#'
#' Lee, E., Chuang, H. Y., Kim, J. W., Ideker, T., & Lee, D. (2008). Inferring
#' pathway activity toward precise disease classification.
#' *PLoS computational biology*, 4(11), e1000217.
#' [doi: 10.1371/journal.pcbi.1000217](https://doi.org/10.1371/journal.pcbi.1000217).
#'
#' @seealso
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
#'
#' @export

compute_zscore <- function(expr_data, signatures) {

    single_sig_zscore <- function(dataset, genes) {
        filt_data <- dataset[genes, ]
        scaled_data <- base::scale(t(filt_data), center = TRUE, scale = TRUE)
        zscore_vec <- base::rowSums(scaled_data) / sqrt(length(genes))
        tibble::enframe(zscore_vec, name = "sample_id", value = "zscore")
    }

    result <- base::lapply(
        signatures,
        FUN = function(genes) {
            single_sig_zscore(dataset = expr_data, genes = genes)
        }
    )
    result <- dplyr::bind_rows(result, .id = "signature_id")
    tidyr::pivot_wider(result,
                       id_cols = sample_id,
                       names_from = signature_id,
                       values_from = zscore)

}
