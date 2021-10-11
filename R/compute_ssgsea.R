#' Single sample GSEA
#'
#' Obtain single sample GSEA scores for a list of gene signatures,
#' as implemented in *Barbie et al. (2009)*. This method is a single sample
#' generalization of the original Gene Set Enrichment Analysis (GSEA) procedure
#' (see *Subramanian et al. (2005)*).
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param signatures A list of named character vectors corresponding to gene signatures.
#' @param norm Type of normalization affecting the single sample scores. Can be one of:
#'  * `raw`, obtain raw scores (_default_);
#'  * `separate`, normalize raw scores across samples for each signature separately.
#'  * `all`, normalize raw scores both across samples and signatures.
#' @param alpha Exponent in the running sum of the score calculation which weights
#'  the gene ranks. Default to \eqn{\alpha = 0.25}. GSEA classic uses \eqn{\alpha = 1}.
#'
#' @return A tibble with the first column indicating sample identifiers (`sample_id`)
#'  and all the other columns ....
#'
#' @references
#' Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F.,
#' Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C., Fröhling, S., Chan, E. M.,
#' Sos, M. L., Michel, K., Mermel, C., Silver, S. J., Weir, B. A., Reiling, J. H.,
#' Sheng, Q., Gupta, P. B., … Hahn, W. C. (2009). Systematic RNA interference
#' reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, 462(7269),
#' 108–112. [doi: 10.1038/nature08460](https://doi.org/10.1038/nature08460).
#'
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation
#' analysis for microarray and RNA-seq data. *BMC bioinformatics*, 14, 7.
#' [doi :10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7).
#'
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#' Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.,
#' & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based
#' approach for interpreting genome-wide expression profiles.
#' *Proceedings of the National Academy of Sciences of the United States of America*,
#' 102(43), 15545–15550. [doi: 10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).
#'
#' @seealso
#'
#' @examples
#' set.seed(123)
#'
#' test_expr <- data.frame(
#'     sample1 = rnorm(10, mean = 8, sd = 2),
#'     sample2 = rnorm(10, mean = 7, sd = 2),
#'     sample3 = rnorm(10, mean = 8, sd = 3),
#'     row.names = paste0("gene", 1:10)
#' )
#'
#' test_list <- list(
#'     gs1 = paste0("gene", c(1, 2, 6)),
#'     gs2 = paste0("gene", 1:5),
#'     gs3 = paste0("gene", c(3, 1, 9)),
#'     gs4 = paste0("gene", c(1, 3, 9))
#' )
#'
#' compute_ssgsea(test_expr, test_list)
#'
#' @export
compute_ssgsea <- function(expr_data, signatures, norm = "raw", alpha = 0.25) {

    es_ssgsea <- function(sample_data, geneset) {
        n_in <- length(geneset)
        n_out <- nrow(expr_data) - n_in
        df <- data.frame(ranks = sample_data)
        df <- df[order(df$ranks, decreasing = TRUE), , drop = FALSE]
        df$prob_in <- 0
        df$prob_out <- 0
        df[rownames(df) %in% geneset, "prob_in"] <- df[rownames(df) %in% geneset, "ranks"]
        df[!rownames(df) %in% geneset, "prob_out"] <- 1
        df["prob_in"] <- cumsum(df["prob_in"]^alpha) / sum(df["prob_in"]^alpha)
        df["prob_out"] <- cumsum(df["prob_out"] / n_out)
        sum(df$prob_in - df$prob_out)
    }

    single_sig_ssgsea <- function(dataset, genes) {
        rank_data <- apply(dataset, MARGIN = 2,
                           FUN = rank, na.last = "keep", ties.method = "average")
        ssgsea_vec <- apply(rank_data, MARGIN = 2, FUN = es_ssgsea, geneset = genes)
        tibble::enframe(ssgsea_vec, name = "sample_id", value = "ssgsea")
    }

    result <- base::lapply(
        signatures,
        FUN = function(genes) {
            single_sig_ssgsea(dataset = expr_data, genes = genes)
        }
    )

    if (norm == "separate") {
        result <- lapply(
            result,
            FUN = function(dataset) {
                dplyr::mutate(
                    dataset,
                    ssgsea = (ssgsea - range(ssgsea)[1]) /
                        (range(ssgsea)[2] - range(ssgsea)[1])
                )
            }
        )
    }

    result <- dplyr::bind_rows(result, .id = "signature_id")

    if (norm == "all") {
        result <- dplyr::mutate(
            result,
            ssgsea = ssgsea / (range(ssgsea)[2] - range(ssgsea)[1])
        )
    }

    tidyr::pivot_wider(result,
                       id_cols = sample_id,
                       names_from = signature_id,
                       values_from = ssgsea)

}
