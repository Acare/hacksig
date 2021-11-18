#' Score samples by gene signatures
#'
#' @description
#' Compute gene signature single sample scores in one of different ways.
#'   You can choose to apply either the *original* procedure or one of three single
#'   sample scoring methods: the combined z-score (*Lee et al., 2008*), the single
#'   sample GSEA (*Barbie et al., 2009*) or the singscore method (*Foroutan et al., 2018*).
#' @details
#' `hack_sig()`
#' For `"original"` method, it is intended the procedure used in the original
#'   publication by the authors for computing the signature score. For example,
#'   it can be a weighted sum of gene expression values and model coefficients or
#'   a simple average of gene expression values for a signature.
#'   `hack_sig()` can compute signature scores with the original method only if
#'   this is a relatively simple procedure (e.g weighted sum or average expression).
#'   For more complex methods
#' @inheritSection ss_methods Algorithm
#' @param signatures It can be a list of signatures or a character vector indicating
#'   keywords for a group of signatures. The default (`"all"`) will cause the
#'   function to compute single sample scores for all the signatures implemented
#'   in `hacksig`.
#' @param method A character string specifying which method to use for computing
#'   the single sample score for each signature. You can choose one of:
#'
#'   * `original`, the original method used by the authors of the signature;
#'   * `zscore`, the combined z-score method implemented in `compute_zscore()`;
#'   * `ssgsea`, the single sample GSEA method implemented in `compute_ssgsea()`;
#'   * `singscore`, the singscore method implemented in `compute_singscore()`;
#' @inherit ss_methods params return references
#' @examples
#' # to obtain raw single sample GSEA scores for all signatures run:
#' hack_sig(test_expr, method = "ssgsea")
#'
#' # to obtain normalized scores, instead run:
#' hack_sig(test_expr, method = "ssgsea", sample_norm = "separate")
#'
#' # You can also change the exponent of the ssGSEA running sum with:
#' hack_sig(test_expr, method = "ssgsea", sample_norm = "separate", alpha = 0.5)
#' @seealso [get_sig_info()] to get valid keywords for signatures.
#'
#'   [check_sig()] for getting information about gene signatures..
#' @export
hack_sig <- function(expr_data, signatures = "all", method = "original",
                     direction = "none", sample_norm = "raw", rank_norm = "none",
                     alpha = 0.25) {
    compute_ss_method <- function(sigs, ss_method) {
        switch (ss_method,
                original = ,
                ssgsea = compute_ssgsea(expr_data, sigs,
                                        sample_norm = sample_norm, rank_norm = rank_norm,
                                        alpha = alpha),
                zscore = compute_zscore(expr_data, sigs),
                singscore = compute_singscore(expr_data, sigs, direction = direction),
                stop("Valid choices for 'method' are 'original', 'zscore', 'ssgsea', 'singscore'",
                     call. = FALSE)
        )
    }
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
    if (is.list(signatures)) {
        signatures <- lapply(signatures, FUN = unique)
        signatures <- signatures[lapply(signatures, length) > 1]
        if (is.null(names(signatures)) == TRUE) {
            names(signatures) <- paste0("sig", seq_along(signatures))
        }
        check_info <- check_sig(expr_data = expr_data, signatures = signatures)
        check_info <- check_info[check_info$n_present == 0,]
        if (nrow(check_info) > 0) {
            signatures <- signatures[!names(signatures) %in% check_info$signature_id]
            rlang::warn(
                rlang::format_error_bullets(
                    c("i" = "No genes are present in 'expr_data' for the following signatures:",
                      stats::setNames(check_info$signature_id, rep_len("x", nrow(check_info))))
                )
            )
        }
        compute_ss_method(signatures, method)
        }
    else if (is.character(signatures)) {
        sig_data <- signatures_data
        signatures <- paste0(signatures, collapse = "|")
        if (signatures != "all") {
            sig_data <- sig_data[grep(signatures, sig_data$signature_keywords,
                                      ignore.case = TRUE), ]
            if (nrow(sig_data) == 0) {
                stop("Provided keyword in 'signatures' does not match any class of signature.",
                     call. = FALSE)
            }
        }
        check_info <- check_sig(expr_data = expr_data, signatures = signatures)
        check_info <- check_info[check_info$n_present == 0, ]
        if (nrow(check_info) > 0) {
            sig_data <- sig_data[!sig_data$signature_id %in% check_info$signature_id, ]
            rlang::warn(
                rlang::format_error_bullets(
                    c("i" = "No genes are present in 'expr_data' for the following signatures:",
                      stats::setNames(check_info$signature_id, rep_len("x", nrow(check_info))))
                )
            )
        }
        sig_list <- lapply(split(sig_data[, c("signature_id", "gene_symbol")],
                                 sig_data$signature_id),
                           FUN = `[[`, 2)
        sig_list <- sig_list[lapply(sig_list, length) > 1]
        if (method != "original") {
            compute_ss_method(sig_list, method)
        } else {
            sig_data <- sig_data[!grepl("cinsarc|estimate|ips", sig_data$signature_keywords), ]
            sig_list <- sig_list[!grepl("cinsarc|estimate|ips", names(sig_list))]
            method_list <- tibble::deframe(
                sig_data[!duplicated(sig_data[, c("signature_id", "signature_method")]),
                         c("signature_id", "signature_method")]
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
                    temp <- compute_ssgsea(expr_data, sig_list[i],
                                           sample_norm = sample_norm, rank_norm = rank_norm,
                                           alpha = alpha)
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
