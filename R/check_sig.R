#' Check signatures feasibility
#'
#' @description
#' `check_sig()` is a helper function that permits to obtain useful information
#' prior to testing your gene signatures.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'   row names and samples as columns.
#' @param signatures It can be a list of signatures or a character string indicating
#'   a keyword for a group of signatures (e.g. "immune" or "ifng" for immune and
#'   interferon gamma signatures respectively). The default (`all`) will cause the
#'   function to check for all the signatures implemented in `hacksig`.
#'
#' @return A tibble with a number of rows equal to the number of input signatures
#'   and four columns: `signature_id`, `n_genes`, `n_present` and `frac_present`.
#'
#'   `n_genes` gives the number of genes composing a signature.
#'   `n_present` and `frac_present` are the number and fraction of genes in a
#'   signature which are present in the expression matrix `expr_data`, respectively.
#'
#' @examples
#' check_sig(test_expr)
#' check_sig(test_expr, hallmark_list)
#'
#' @importFrom rlang .data
#' @export
check_sig <- function(expr_data, signatures = "all") {
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
    if (is.list(signatures)) {
        signatures <- lapply(signatures, FUN = unique)
        if (is.null(names(signatures)) == TRUE) {
            names(signatures) <- paste0("sig", seq_along(signatures))
        }
        sig_data <- tidyr::unnest(
            tibble::enframe(signatures, name = "signature_id", value = "gene_symbol"),
            cols = "gene_symbol"
        )
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
    }
    sig_data_group <- dplyr::group_by(sig_data[, c("signature_id", "gene_symbol")],
                                      .data$signature_id)
    sig_data_group <- dplyr::mutate(
        sig_data_group,
        n_genes = length(.data$gene_symbol),
        n_present = length(intersect(.data$gene_symbol, rownames(expr_data))),
        frac_present = .data$n_present / .data$n_genes
    )
    result <- dplyr::ungroup(
        dplyr::distinct(
            sig_data_group[, c("signature_id", "n_genes", "n_present", "frac_present")]
        )
    )
    dplyr::arrange(result, -.data$frac_present)
}
