#' Score samples by gene signatures
#'
#' `hack_sig()` is the main function of the package, which computes single sample
#'  scores in one of different ways. You can choose to apply the default method
#'  for a signature or you can choose one of three single sample scoring methods:
#'  **z-score** (*Lee et al., 2008*), **ssGSEA** (*Barbie et al., 2009*) or
#'  **singscore** (*Foroutan et al., 2018*).
#'
#' @param method A character string specifying which method to use for computing
#'  the single sample score for each signature. You can choose one of:
#'
#'  * `default`, the original method used by the authors of the signature;
#'  * `zscore`, the z-score method implemented in `compute_zscore()`;
#'  * `ssgsea`, the single sample GSEA method implemented in `compute_ssgsea()`;
#'  * `singscore`, the singscore method implemented in `compute_singscore()`;
#' @inherit ss_methods params return references
#'
#' @examples
#' # to obtain raw single sample GSEA scores for all signatures run:
#' hack_sig(test_expr, method = "ssgsea")
#' # to obtain normalized scores, instead run:
#' hack_sig(test_expr, method = "ssgsea", norm = "separate")
#' # You can also change the exponent of the ssGSEA running sum with:
#' hack_sig(test_expr, method = "ssgsea", norm = "separate", alpha = 0.5)
#'
#' @export
hack_sig <- function(expr_data, signatures = "all", method = "default",
                     direction = "none", norm = "raw", alpha = 0.25) {
    compute_ss_method <- function(ss_method) {
        switch (ss_method,
                default = ,
                ssgsea = compute_ssgsea(expr_data, signatures, norm = norm, alpha = alpha),
                zscore = compute_zscore(expr_data, signatures),
                singscore = compute_singscore(expr_data, signatures, direction = direction),
                stop("Valid choices for 'method' are 'default', 'zscore', 'ssgsea', 'singscore'",
                     call. = FALSE)
        )
    }
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
    if (is.list(signatures)) {
        signatures <- lapply(signatures, FUN = unique)
        if (is.null(names(signatures)) == TRUE) {
            names(signatures) <- paste0("sig", seq_along(signatures))
        }
        compute_ss_method(method)
        }
    else if (is.character(signatures)) {
        sig_data <- hacksig::signatures_data
        if (signatures != "all") {
            sig_data <- sig_data[grep(signatures, sig_data$signature_id), ]
            if (nrow(sig_data) == 0) {
                stop("Provided keyword in 'signatures' does not match any class of signature.",
                     call. = FALSE)
            }
        }
        sig_list <- lapply(split(sig_data[, c("signature_id", "gene_symbol")],
                                 sig_data$signature_id),
                           FUN = `[[`, 2)
        if (method != "default") {
            compute_ss_method(method)
        } else {
            method_list <- tibble::deframe(
                sig_data[!duplicated(sig_data[, c("signature_id", "method")]),
                         c("signature_id", "method")]
            )
            method_list <- method_list[match(names(sig_list), names(method_list))]
            weight_list <- lapply(split(sig_data[, c("signature_id", "gene_weight")],
                                        sig_data$signature_id),
                                  FUN = `[[`, 2)
            result <- vector("list", length = length(sig_list))
            for (i in names(method_list)) {
                if (grepl("weighted_sum", method_list[[i]])) {
                    temp <- tibble::enframe(
                        colSums(expr_data[sig_list[[i]], ] * weight_list[[i]],
                                na.rm = TRUE),
                        name = "sample_id",
                        value = "sig_score"
                    )
                    temp$signature_id <- i
                } else if (grepl("mean", method_list[[i]])) {
                    temp <- tibble::enframe(
                        colMeans(expr_data[sig_list[[i]], ], na.rm = TRUE),
                        name = "sample_id",
                        value = "sig_score"
                    )
                    temp$signature_id <- i
                } else {
                    temp <- compute_ssgsea(expr_data, sig_list[i], norm = norm, alpha = alpha)
                    temp$sig_score <- temp[, i, drop = TRUE]
                    temp[, i] <- NULL
                    temp$signature_id <- i
                }
                result[[i]] <- temp
            }
            result <- dplyr::bind_rows(result)

            tidyr::pivot_wider(result,
                               id_cols = "sample_id",
                               names_from = "signature_id",
                               values_from = "sig_score")
        }
    } else stop("Argument 'signatures' must be a named list of gene signatures or a string with a keyword.",
                call. = FALSE)
}
