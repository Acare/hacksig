#' Check signatures feasibility
#'
#' @description
#' `check_sig()` is a helper function that shows useful information about signatures
#'   that you want to test on your gene expression matrix.
#' @param signatures It can be a list of signatures or a character vector indicating
#'   keywords for a group of signatures. The default (`"all"`) will cause the
#'   function to check for all the signatures implemented in `hacksig`.
#' @inheritParams hack_estimate
#' @return A tibble with a number of rows equal to the number of input signatures
#'   and five columns:
#'
#'   * `n_genes` gives the number of genes composing a signature;
#'   * `n_present` and `frac_present` are the number and fraction of genes in a
#'     signature which are present in `expr_data`, respectively;
#'   * `missing_genes` returns a named list of missing gene symbols for each signature.
#' @examples
#' check_sig(test_expr)
#' check_sig(test_expr, "estimate")
#' @importFrom rlang .data
#' @seealso [get_sig_info()], [hack_sig()]
#' @export
check_sig <- function(expr_data, signatures = "all") {
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
    if (is.list(signatures) == TRUE) {
        signatures <- lapply(signatures, FUN = unique)
        if (is.null(names(signatures)) == TRUE) {
            names(signatures) <- paste0("sig", seq_along(signatures))
        }
        sig_data <- tidyr::unnest(
            tibble::enframe(signatures, name = "signature_id", value = "gene_symbol"),
            cols = "gene_symbol"
        )
    }
    else if (is.character(signatures) == TRUE) {
        sig_data <- signatures_data
        signatures <- paste0(signatures, collapse = "|")
        if (signatures != "all") {
            sig_data <- sig_data[grep(signatures, sig_data$signature_keywords), ]
            if (nrow(sig_data) == 0) {
                stop("Provided keyword in 'signatures' does not match any class of signature.",
                     call. = FALSE)
            }
        }
    }
    sig_data_group <- dplyr::group_by(sig_data[, c("signature_id", "gene_symbol")],
                                      .data$signature_id)
    sig_data_group <- dplyr::mutate(
        sig_data_group,
        n_genes = length(.data$gene_symbol),
        n_present = length(intersect(.data$gene_symbol, rownames(expr_data))),
        frac_present = .data$n_present / .data$n_genes,
        missing_genes = list(setdiff(.data$gene_symbol, rownames(expr_data))),
        missing_genes = stats::setNames(.data$missing_genes, .data$signature_id)
    )
    keep_cols <- c("signature_id", "n_genes", "n_present", "frac_present", "missing_genes")
    result <- unique(
        dplyr::ungroup(
            sig_data_group[, keep_cols]
        )
    )
    dplyr::arrange(result, -.data$frac_present)
}
