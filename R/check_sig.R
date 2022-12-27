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
#'   * `signature_id`, the signature ID;
#'   * `n_genes`, the number of genes composing a signature;
#'   * `n_present` and `frac_present`, the number and fraction of genes in a
#'     signature which are present in `expr_data`, respectively;
#'   * `missing_genes`, the missing gene symbols for each signature.
#' @examples
#' check_sig(test_expr)
#' check_sig(test_expr, "estimate")
#' @seealso [get_sig_info()], [hack_sig()]
#' @importFrom data.table .N `:=`
#' @export
check_sig <- function(expr_data, signatures = "all") {
    if (is.matrix(expr_data)) {
        expr_data <- as.data.frame(expr_data)
    }
    if (is.list(signatures)) {
        signatures <- lapply(signatures, FUN = unique)
        if (is.null(names(signatures))) {
            names(signatures) <- paste0("sig", seq_along(signatures))
        }
        sig_data <- data.table::rbindlist(
            tibble::enframe(signatures, name = "signature_id", value = "gene_symbol"),
            idcol = "gene_symbol"
        )
    }
    else if (is.character(signatures)) {
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
    data.table::setDT(sig_data)
    result <- sig_data[
        ,
        .(n_genes = .N,
          n_present = length(intersect(gene_symbol, rownames(expr_data))),
          missing_genes = list(setdiff(gene_symbol, rownames(expr_data)))),
        by = "signature_id"
    ][
        ,
        frac_present := n_present / n_genes
    ]
    data.table::setcolorder(result,
                            c("signature_id", "n_genes", "n_present", "frac_present", "missing_genes"))
    data.table::setorderv(result, cols = c("frac_present", "n_genes"), order = -1)
    tibble::as_tibble(result)
}
