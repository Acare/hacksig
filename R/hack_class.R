#' Hack signature classes
#'
#' @description
#' `hack_class()` is supposed to be used in combination after [hack_sig()] in
#'   order to classify your samples in one of two or more signature classes.
#' @param sig_data A tibble result of a call to [hack_sig()].
#' @param cutoff A character specifying which function to use to categorize
#'   samples by signature scores. Can be one of:
#'
#'   * `"original"` (default), apply the original publication method; if
#'     categorization is not expected, the median score is used as a threshold;
#'   * `"mean"/"median"`, samples will be classified as `"low"` or `"high"` with
#'     respect to the mean/median signature score, respectively;
#'   * `"tertiles"`, samples will be classified as `"<= T1"` (score lower than
#'     first tertile), `"(T1, T2]"` (score between first and second tertiles),
#'     `"> T2"` (score higher than second tertile);
#'   * `"quartiles"`, samples will be classified as `"<= Q1"` (score lower than
#'     first quartile), `"(Q1, Q2]"` (score between first and second quartiles),
#'     `"(Q2, Q3]"` (score between second and third quartiles), `"> Q3"` (score
#'     higher than third quartile).
#' @return A tibble with the same dimension as `sig_data`, a column `sample_id`
#'   indicating sample identifiers and one column for each input signature giving
#'   signature classes.
#' @examples
#' library(dplyr)
#' hack_sig(test_expr, "immune") %>% hack_class()
#' @seealso [hack_sig()]
#' @export
hack_class <- function(sig_data, cutoff = "original") {
    sig_info <- signatures_data
    sig_info <- sig_info[sig_info$signature_id %in% names(sig_data), ]
    if (cutoff == "original") {
        method_list <- tibble::deframe(
            sig_info[!duplicated(sig_info[, c("signature_id", "signature_method")]),
                     c("signature_id", "signature_method")]
        )
        method_list <- gsub(".*\\|", "", method_list)
        result <- sig_data
        for (i in names(method_list)) {
            if (grepl("mean", method_list[[i]])) {
                result[i] <- ifelse(result[i] < mean(result[[i]], na.rm = TRUE),
                                      "low", "high")
            } else if (grepl("\\d", method_list[[i]])) {
                result[i] <- ifelse(result[i] < as.numeric(method_list[[i]]),
                                      "low", "high")
            } else {
                result[i] <- ifelse(result[i] < stats::median(result[[i]], na.rm = TRUE),
                                      "low", "high")
            }
        }
        result
    } else if (cutoff %in% c("mean", "median")) {
        result <- apply(
            tibble::column_to_rownames(sig_data, "sample_id"),
            MARGIN = 2,
            FUN = function(x) ifelse(x < do.call(cutoff, list(x = x, na.rm = TRUE)),
                                     "low", "high")
        )
        tibble::as_tibble(result, rownames = "sample_id")
    } else if (cutoff == "tertiles") {
        dplyr::mutate(
            sig_data,
            dplyr::across(
                .cols = -"sample_id",
                .fns = function(x) {
                    cut(x,
                        breaks = stats::quantile(x, probs = seq(from = 0, to = 1, by = 1 / 3)),
                        labels = c("<= T1", "(T1, T2]", "> T2"),
                        include.lowest = TRUE)
                }
            )
        )
    } else if (cutoff == "quartiles") {
        dplyr::mutate(
            sig_data,
            dplyr::across(
                .cols = -"sample_id",
                .fns = function(x) {
                    cut(x,
                        breaks = stats::quantile(x, probs = seq(from = 0, to = 1, by = 0.25)),
                        labels = c("<= Q1", "(Q1, Q2]", "(Q2, Q3]", "> Q3"),
                        include.lowest = TRUE)
                }
            )
        )
    } else stop("Valid choices for 'cutoff' are 'original', 'mean', 'median', 'tertiles' and 'quartiles'",
                call. = FALSE)
}
