#' Display available gene signatures
#'
#' @description
#' `get_sig_info()` returns information about all gene signatures implemented in
#'   `hacksig`.
#' @return A tibble with one row per signature and four columns:
#'
#'   * `signature_id`, a unique identifier associated to a signature;
#'   * `signature_keywords`, valid keywords to use in the `signatures` argument of
#'   `hack_sig()` and `check_sig()` as well as in the `keywords` argument of
#'   `get_sig_genes()`;
#'   * `publication_doi`, the original publication DOI;
#'   * `description`, a brief description about the signature.
#' @seealso [check_sig()], [hack_sig()], [get_sig_genes()]
#' @examples
#' get_sig_info()
#' @importFrom data.table `:=`
#' @export
get_sig_info <- function() {
    signature_keywords = NULL # due to NSE notes in R CMD check
    sig_data <- unique(
        signatures_data[, c("signature_id", "signature_keywords", "publication_doi", "description")]
    )
    data.table::setDT(sig_data)
    sig_data[, signature_keywords := paste0(signature_keywords, collapse = "|"), by = "signature_id"]
    sig_data <- unique(sig_data)
    sig_data[
        ,
        signature_keywords := paste0(sort(unique(strsplit(signature_keywords, "\\|")[[1]])), collapse = "|"),
        by = "signature_id"
    ]
    data.table::setorderv(sig_data, cols = "signature_id", order = 1)
    tibble::as_tibble(sig_data)
}
