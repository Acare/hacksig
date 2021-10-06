#' ESTIMATE
#'
#' Obtain Immune, Stroma, ESTIMATE and Purity scores from a cohort of samples.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#'
#' @return A data frame with
#' @export
#'
#' @examples
#' m <- matrix()
hack_estimate <- function(expr_data) {
    # calcola ssgsea grezzi per immune e stroma
    # calcola estimate facendo la somma
    # calcola purity con formula e mettendo NA quelli < 0
    # ritorna dataframe con campioni, immune, stroma, estimate e purity scores
}
