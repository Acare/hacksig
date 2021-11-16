#' Get keywords associated to signatures
#'
#' @description
#' Obtain valid keywords for the argument `signatures` of `hack_sig()` and
#' `check_sig()` in order to filter for specific gene signatures.
#' Some of the returned strings are redundant, such as `"interferon gamma"` and `"ifng"`.
#' @return A character vector with valid keywords.
#' @seealso [check_sig()], [hack_sig()]
#' @examples
#' get_sig_keywords()
#' @export
get_sig_keywords <- function() {
    sig_data <- hacksig::signatures_data
    sort(
        unique(
            unlist(
                strsplit(sig_data$signature_keyword, split = "|", fixed = TRUE)
            )
        )
    )
}
