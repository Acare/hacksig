#' Display available gene signatures
#'
#' @description
#' `get_sig_info()` returns information about all gene signatures implemented in
#'   `hacksig`.
#' @return A tibble with one row per signature and four columns:
#'
#'   * `signature_id`;
#'   * `signature_keywords`, valid keywords to use in the `signatures` argument of
#'   `hack_sig()` and `check_sig()`;
#'   * `publication_doi`, the original publication DOI;
#'   * `description`, a brief description about the signature usage.
#' @seealso [check_sig()], [hack_sig()]
#' @examples
#' get_sig_info()
#' @importFrom rlang .data
#' @export
get_sig_info <- function() {
    sig_data <- signatures_data
    keep_cols <- c("signature_id", "signature_keywords", "publication_doi", "description")
    sig_data <- unique(sig_data[, keep_cols])
    dplyr::arrange(sig_data, .data$signature_id)
}
