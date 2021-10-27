#' Hack the ESTIMATE scores
#'
#' @description
#' Obtain *Immune*, *Stroma*, *ESTIMATE* and *Purity* scores from a cohort of samples,
#' as implemented in *Yoshihara et al., 2013*.
#'
#' @details
#' ESTIMATE (*Estimation of STromal and Immune cells in MAlignant Tumours using
#' Expression data*) is a method that, through the use of gene expression signatures,
#' infers the fraction of stromal and immune cells in tumor samples. Stromal and
#' immune cells play an important biological role, but can lead a perturbation in
#' gene expression signals. The ESTIMATE score is able to predict the tumour
#' purity in tumour tissue. It was developed on The Cancer Genome Atlas data
#' across 11 different tumor types (Bladder urothelial carcinoma, Breast cancer,
#' Colon and rectal adenocarcinoma, Glioblastoma multiforme, Head and neck
#' squamous cell carcinoma, Clear cell renal cell carcinoma, Lung adenocarcinoma,
#' Lung squamous cell carcinoma, Ovarian serous cystadenocarcinoma, Uterine corpus
#' endometrial carcinoma), and confirmed using 3,809 transcriptional profiles
#' present in available public datasets.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#'
#' @return A tibble with one row for each sample in `expr_data` and five columns:
#'   `sample_id`, `immune_score`, `stroma_score`, `estimate_score` and `purity_score`.
#'
#' @references
#' Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F.,
#' Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C., Fröhling, S., Chan, E. M.,
#' Sos, M. L., Michel, K., Mermel, C., Silver, S. J., Weir, B. A., Reiling, J. H.,
#' Sheng, Q., Gupta, P. B., … Hahn, W. C. (2009). Systematic RNA interference
#' reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, 462(7269),
#' 108–112. [doi: 10.1038/nature08460](https://doi.org/10.1038/nature08460).
#'
#' Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R.,
#' Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A.,
#' Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013).
#' Inferring tumour purity and stromal and immune cell admixture from
#' expression data. *Nature communications*, 4, 2612.
#' [doi: 10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612)
#'
#' @examples
#' hack_estimate(test_expr)
#' @export
hack_estimate <- function(expr_data) {
    sig_data <- hacksig::signatures_data
    estimate_df <- sig_data[sig_data$signature_id == "estimate",
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
