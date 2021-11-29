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
#' @importFrom rlang .data
#' @export
get_sig_info <- function() {
    sig_data <- signatures_data
    keep_cols <- c("signature_id", "signature_keywords", "publication_doi", "description")
    sig_data <- unique(sig_data[, keep_cols])
    split_kwrds <- dplyr::mutate(
        dplyr::rowwise(sig_data),
        signature_keywords = list(strsplit(.data$signature_keywords, "\\|"))
    )
    split_kwrds <- dplyr::summarise(
        dplyr::group_by(split_kwrds, .data$signature_id),
        signature_keywords = list(unique(unlist(lapply(.data$signature_keywords, `[[`, 1))))
    )
    split_kwrds <- dplyr::mutate(
        dplyr::group_by(split_kwrds, .data$signature_id),
        signature_keywords = paste0(sort(unlist(.data$signature_keywords)),
                                    collapse = "|")
    )
    sig_data <- dplyr::left_join(split_kwrds,
                                 unique(sig_data[, keep_cols[-2]]),
                                 by = "signature_id")
    dplyr::arrange(dplyr::ungroup(sig_data), .data$signature_id)
}
