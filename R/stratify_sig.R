#' Stratify samples into classes
#'
#' @description
#' `stratify_sig()` is supposed to be used in combination after [hack_sig()] in
#'   order to classify your samples in one of two or more signature classes.
#' @param sig_data A tibble result of a call to [hack_sig()].
#' @param cutoff A character specifying which function to use to categorize
#'   samples by signature scores. Can be one of:
#'
#'   * `"original"` (default), apply the original publication method; if
#'     categorization is not expected, the median score is used as a threshold;
#'   * `"mean"`/`"median"`, samples will be classified as `"low"` or `"high"` with
#'     respect to the mean/median signature score, respectively;
#'   * `"quantile"`, samples will be classified into signature score quantiles;
#' @param probs A numeric vector of probabilities with values in $[0, 1]$ to use in
#'   combination with `cutoff = "quantile"`. By default, it correspond to quartiles
#'   (`c(0, 0.25, 0.5, 0.75, 1)`).
#' @return A tibble with the same dimension as `sig_data`, having a column `sample_id`
#'   with sample identifiers and one column for each input signature giving
#'   sample classes.
#' @examples
#' scores <- hack_sig(test_expr, "immune")
#' stratify_sig(scores)
#' @seealso [hack_sig()], [stats::quantile()]
#' @export
stratify_sig <- function(sig_data, cutoff = "original", probs = seq(0, 1, 0.25)) {
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
    } else if (cutoff == "quantile") {
        data.table::setDT(sig_data)
        sig_data <- sig_data[
            ,
            (setdiff(colnames(sig_data), "sample_id")) :=
                lapply(.SD,
                       function(x) {
                           paste0(
                               "Q",
                               cut(x,
                                   breaks = stats::quantile(x, probs = probs),
                                   labels = FALSE,
                                   include.lowest = TRUE)
                           )
                       }),
            .SDcols = setdiff(colnames(sig_data), "sample_id")
        ]
        tibble::as_tibble(sig_data)
    } else stop("Valid choices for 'cutoff' are 'original', 'mean', 'median' and 'quantile'",
                call. = FALSE)
}
