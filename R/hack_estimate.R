#' Hack the ESTIMATE scores
#'
#' @description
#' Obtain *Immune*, *Stroma*, *ESTIMATE* and *Tumor Purity* scores from a cohort
#'   of samples, using the method implemented in *Yoshihara et al., 2013*.
#' @details
#' The ESTIMATE (*Estimation of STromal and Immune cells in MAlignant Tumors
#'   using Expression data*) method was developed with the aim to estimate the
#'   fraction of tumor cells in a sample by using gene expression instead of copy
#'   number data. The fundamental assumption of this method is that the tumor
#'   microenvironment is a very rich and dynamic ecosystem, in which immune
#'   infiltrating cells and stroma play a major role. The ESTIMATE score is
#'   defined as the combination (i.e. sum) of immune and stroma scores and can be
#'   thought of as a *"non-tumor score"*. Consequently, a high ESTIMATE enrichment
#'   gives a low tumor purity score and viceversa.
#' @section Algorithm:
#' Raw immune and stromal signatures scores are computed using single sample GSEA
#'   with rank normalization (*Barbie et al., 2009*).
#'   Then, the ESTIMATE score is computed by summing the immune and stroma scores.
#'   Finally, the tumor purity score is obtained with the following formula:
#'   \deqn{Purity = cos(0.6049872018 + 0.0001467884 * ESTIMATE)}
#' @param expr_data A normalized gene expression matrix (or data frame) with
#'   gene symbols as row names and samples as columns.
#' @return A tibble with one row for each sample in `expr_data` and five columns:
#'   `sample_id`, `immune_score`, `stroma_score`, `estimate_score` and `purity_score`.
#' @source [bioinformatics.mdanderson.org/public-software/estimate/](https://bioinformatics.mdanderson.org/public-software/estimate/)
#' @references
#' Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F.,
#' Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C., Fröhling, S., Chan, E. M.,
#' Sos, M. L., Michel, K., Mermel, C., Silver, S. J., Weir, B. A., Reiling, J. H.,
#' Sheng, Q., Gupta, P. B., … Hahn, W. C. (2009). Systematic RNA interference
#' reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, 462(7269),
#' 108–112. \doi{10.1038/nature08460}.
#'
#' Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R.,
#' Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A.,
#' Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013).
#' Inferring tumour purity and stromal and immune cell admixture from
#' expression data. *Nature communications*, 4, 2612. \doi{10.1038/ncomms3612}.
#' @examples
#' hack_estimate(test_expr)
#' @export
hack_estimate <- function(expr_data) {
    sig_data <- signatures_data
    estimate_sigs <- list(
        immune_score = sig_data[sig_data$signature_id == "estimate_immune",
                                "gene_symbol",
                                drop = TRUE],
        stroma_score = sig_data[sig_data$signature_id == "estimate_stromal",
                                "gene_symbol",
                                drop = TRUE]
    )
    result <- compute_ssgsea(expr_data = expr_data, signatures = estimate_sigs,
                             sample_norm = "raw", rank_norm = "rank", alpha = 0.25)
    result$estimate_score <- result$immune_score + result$stroma_score
    result$purity_score <- cos(0.6049872018 + 0.0001467884 * result$estimate_score)
    result[result$purity_score < 0, "purity_score"] <- NA
    result
}
