#' Get signature gene identifiers
#'
#' @description
#' Obtain gene signatures implemented in `hacksig` as a named list.
#'
#' @param keyword A character string indicating a keyword for a group of signatures.
#'   The default (`"all"`) will cause the function to check for all the signatures
#'   implemented in `hacksig`.
#' @return A named list of gene signatures.
#' @seealso [get_sig_keywords()] to get valid keywords for signatures.
#' @examples
#' get_sig_genes()
#' get_sig_genes("estimate")
#' @export
get_sig_genes <- function(keyword = "all") {
    sig_data <- hacksig::signatures_data
    if (keyword != "all") {
        sig_data <- sig_data[grep(keyword, sig_data$signature_keyword,
                                  ignore.case = TRUE), ]
        if (nrow(sig_data) == 0) {
            stop("Provided keyword does not match any class of signature.",
                 call. = FALSE)
        }
    }
    sig_list <- lapply(split(sig_data[, c("signature_id", "gene_symbol")],
                             sig_data$signature_id),
                       FUN = `[[`, 2)
    sig_list <- sig_list[lapply(sig_list, length) > 1]
    lapply(sig_list, sort)
}
