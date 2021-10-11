#' ESTIMATE
#'
#' Obtain Immune, Stroma, ESTIMATE and Purity scores from a cohort of samples.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#'
#' @return A data frame with
#'
#' @references Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R.,
#' Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A.,
#' Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013).
#' Inferring tumour purity and stromal and immune cell admixture from
#' expression data. *Nature communications*, 4, 2612.
#' [doi: 10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612)
#'
#' @examples
#'
#' @export
hack_estimate <- function(expr_data) {
    estimate_df <- signatures_data
    estimate_df <- estimate_df[estimate_df$signature_id == "estimate",
                               c("gene_type", "gene_symbol")]
    estimate_sigs <- list(
        immune_score = estimate_df[estimate_df$gene_type == "immune",
                                   "gene_symbol",
                                   drop = TRUE],
        stroma_score = estimate_df[estimate_df$gene_type == "stromal",
                                   "gene_symbol",
                                   drop = TRUE]
    )
    result <- compute_ssgsea(expr_data = expr_data, signatures = estimate_sigs)
    result$estimate_score <- result$immune_score + result$stroma_score
    result$purity_score <- cos(0.6049872018 + 0.0001467884 * result$estimate_score)
    result[result$purity_score < 0, "purity_score"] <- NA
    result
}
