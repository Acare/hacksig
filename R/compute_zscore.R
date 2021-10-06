#' Single Sample Z-score
#'
#' Obtain single sample z-scores, as implemented in XXX, from a cohort of samples.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param signature A list of character vectors corresponding to gene signatures.
#'
#' @return A data frame with
#' @export
#'
#' @examples
#' m <- matrix()
compute_zscore <- function(expr_data, signatures) {

    single_sig_zscore <- function(dataset, genes) {
        filt_data <- dataset[genes, ]
        scaled_data <- scale(t(filt_data), center = TRUE, scale = TRUE)
        zscore_vec <- rowSums(scaled_data) / sqrt(length(genes))
        zscore_df <-  tibble::enframe(zscore_vec,
                                      name = "sample_id", value = "zscore")
    }

    result <- lapply(signatures, single_sig_zscore,
                     dataset = expr_data)

    result <- bind_rows(result, .id = "signature_id")
    result <- tidyr::pivot_wider(result, sample_id,
                                 names_from = signature_id,
                                 values_from = zscore)
    result

    # filt_data <- expr_data[signatures, ]
    # scaled_data <- scale(t(filt_data), center = TRUE, scale = TRUE)
    # zscore_vec <- rowSums(scaled_data) / sqrt(length(signatures))
    # result <-  tibble::enframe(zscore_vec, name = "sample_id", value = "zscore")
    #
    # result
}
