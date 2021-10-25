#' Score samples by gene signatures
#'
#' `hack_sig()` is the main function of the package, which computes single sample
#'  scores in one of different ways. You can choose to apply the default method
#'  for a signature or you can choose one of three single sample scoring methods:
#'   * __z-score__ (_Lee et al., 2008_);
#'   * __ssGSEA__ (_Barbie et al., 2009_);
#'   * __singscore__ (_Foroutan et al., 2018_).
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param signatures It can be a list of signatures or a character string indicating
#'  a keyword for a group of signatures (e.g. "immune" or "ifng" for immune and
#'   interferon gamma signatures respectively). The default (`all`) will cause the
#'   function to compute single sample scores for all the signatures implemented
#'   in `hacksig`.
#' @param method A character string specifying which method to use for computing
#'  the single sample score for each signature. You can choose one of:
#'
#'  * `default`, the original method used by the authors of the signature;
#'  * `zscore`, the z-score method implemented in `compute_zscore()`;
#'  * `ssgsea`, the single sample GSEA method implemented in `compute_ssgsea()`;
#'  * `singscore`, the singscore method implemented in `compute_singscore()`;
#' @param ... Additional arguments passed to `compute_zscore()`, `compute_ssgsea()`
#'  or `compute_singscore()`, depending on the choice of `method`.
#' @param direction A character specifying the score computation depending on the
#'  direction of the signatures. Can be on of:
#'  * `none`, undirected signatures, that is you don't know whether the genes are
#'   up- or down-regulated (default);
#'  * `up`, all genes in the signature are supposed to be up-regulated;
#'  * `down`, all genes in the signature are supposed to be down-regulated;
#'  * `both`, a signature is composed of up- and down sub-signatures.
#'   You must supply it as a nested list.
#' @param norm A character string specifying the type of normalization affecting
#'  the single sample GSEA scores. Can be one of:
#'   * `raw`, obtain raw scores (_default_);
#'   * `separate`, normalize raw scores across samples for each signature separately.
#'   * `all`, normalize raw scores both across samples and signatures.
#' @param alpha A numeric scalar. Exponent in the running sum of the single sample
#'  GSEA score calculation which weights the gene ranks. Defaults to \eqn{\alpha = 0.25}.
#' @return A tibble with the first column indicating sample identifiers (`sample_id`)
#'  and all the other columns ....
#' @references
#' Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E., Dunn, I. F.,
#' Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C., Fröhling, S., Chan, E. M.,
#' Sos, M. L., Michel, K., Mermel, C., Silver, S. J., Weir, B. A., Reiling, J. H.,
#' Sheng, Q., Gupta, P. B., … Hahn, W. C. (2009). Systematic RNA interference
#' reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, 462(7269),
#' 108–112. [doi: 10.1038/nature08460](https://doi.org/10.1038/nature08460).
#'
#' Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., &
#' Davis, M. J. (2018). Single sample scoring of molecular phenotypes.
#' *BMC bioinformatics*, 19(1), 404.
#' [doi: 10.1186/s12859-018-2435-4](https://doi.org/10.1186/s12859-018-2435-4).
#'
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation
#' analysis for microarray and RNA-seq data. *BMC bioinformatics*, 14, 7.
#' [doi: 10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7).
#'
#' Lee, E., Chuang, H. Y., Kim, J. W., Ideker, T., & Lee, D. (2008). Inferring
#' pathway activity toward precise disease classification.
#' *PLoS computational biology*, 4(11), e1000217.
#' [doi: 10.1371/journal.pcbi.1000217](https://doi.org/10.1371/journal.pcbi.1000217).
#'
#' Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
#' Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S.,
#' & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based
#' approach for interpreting genome-wide expression profiles.
#' *Proceedings of the National Academy of Sciences of the United States of America*,
#' 102(43), 15545–15550. [doi: 10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).
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
#' hack_sig(test_expr, test_list)
#' @export hack_sig
hack_sig <- function(expr_data, signatures = "all", method = "default", ...) {
    compute_ss_method <- function(ss_method) {
        switch (ss_method,
                default = ,
                ssgsea = .compute_ssgsea(expr_data, signatures, ...),
                zscore = .compute_zscore(expr_data, signatures),
                singscore = .compute_singscore(expr_data, signatures, ...),
                stop("Valid choices for 'method' are 'default', 'zscore', 'ssgsea', 'singscore'",
                     call. = FALSE)
        )
    }
    if (is.list(signatures)) {
        # check if unique gene identifiers
        # check if named list, otherwise give names (sig1, sig2, etc...)
        compute_ss_method(method)
        }
    else if (is.character(signatures)) {
        sig_data <- hacksig::signatures_data
        if (signatures != "all") {
            sig_data <- sig_data[grep(signatures, signatures_data$signature_id), ]
            if (nrow(sig_data) == 0) {
                stop("Provided keyword in 'signatures' does not match any class of signature.",
                     call. = FALSE)
            }
        }
        sig_list <- lapply(split(sig_data[, c("signature_id", "gene_symbol")],
                                 sig_data$signature_id),
                           FUN = `[[`, 2)
        if (method != "default") {
            compute_ss_method(method)
        } else {
            method_list <- tibble::deframe(
                sig_data[!duplicated(sig_data[, c("signature_id", "method")]),
                         c("signature_id", "method")]
            )
            method_list <- method_list[match(names(sig_list), names(method_list))]
            weight_list <- lapply(split(sig_data[, c("signature_id", "gene_weight")],
                                        sig_data$signature_id),
                                  FUN = `[[`, 2)
            result <- vector("list", length = length(sig_list))
            for (i in names(method_list)) {
                if (grepl("weighted_sum", method_list[[i]])) {
                    temp <- tibble::enframe(
                        colSums(expr_data[sig_list[[i]], ] * weight_list[[i]],
                                na.rm = TRUE),
                        name = "sample_id",
                        value = "sig_score"
                    )
                    temp$signature_id <- i
                } else if (grepl("mean", method_list[[i]])) {
                    temp <- tibble::enframe(
                        colMeans(expr_data[sig_list[[i]], ], na.rm = TRUE),
                        name = "sample_id",
                        value = "sig_score"
                    )
                    temp$signature_id <- i
                } else {
                    temp <- .compute_ssgsea(expr_data, sig_list[i], ...)
                    temp$sig_score <- temp[, i, drop = TRUE]
                    temp[, i] <- NULL
                    temp$signature_id <- i
                }
                result[[i]] <- temp
            }
            result <- dplyr::bind_rows(result)

            tidyr::pivot_wider(result,
                               id_cols = "sample_id",
                               names_from = "signature_id",
                               values_from = "sig_score")
        }
    } else stop("Argument 'signatures' must be a named list of gene signatures or a string with a keyword.",
                call. = FALSE)
}

.compute_zscore <- function(expr_data, signatures) {
    single_sig_zscore <- function(dataset, genes) {
        filt_data <- dataset[rownames(dataset) %in% genes, ]
        scaled_data <- scale(t(filt_data), center = TRUE, scale = TRUE)
        zscore_vec <- rowSums(scaled_data, na.rm = TRUE) / sqrt(ncol(scaled_data))
        tibble::enframe(zscore_vec, name = "sample_id", value = "zscore")
    }

    result <- future.apply::future_lapply(
        X = signatures,
        FUN = function(genes) {
            single_sig_zscore(dataset = expr_data, genes = genes)
        }
    )
    result <- dplyr::bind_rows(result, .id = "signature_id")

    tidyr::pivot_wider(result,
                       id_cols = "sample_id",
                       names_from = "signature_id",
                       values_from = "zscore")

}

.compute_ssgsea <- function(expr_data, signatures, norm = "raw", alpha = 0.25) {
    n_genes <- nrow(expr_data)

    es_ssgsea <- function(sample_data, geneset) {
        n_in <- length(intersect(geneset, rownames(expr_data)))
        n_out <- n_genes - n_in
        # check rank options as in GSVA::gsva()
        # rank_vec <- rank(sample_data, na.last = "keep", ties.method = "average")
        rank_vec <- rank(sample_data)
        rank_vec <- rank_vec[order(rank_vec, decreasing = TRUE)]
        prob_in <- rep.int(0, length(rank_vec))
        prob_out <- rep.int(0, length(rank_vec))
        is_geneset_vec <- names(rank_vec) %in% geneset
        prob_in[is_geneset_vec] <- rank_vec[is_geneset_vec]
        prob_out[!is_geneset_vec] <- 1
        prob_in <- cumsum(prob_in^alpha) / sum(prob_in^alpha)
        prob_out <- cumsum(prob_out / n_out)
        sum(prob_in - prob_out)
    }

    single_sig_ssgsea <- function(dataset, genes) {
        ssgsea_vec <- apply(dataset, MARGIN = 2,
                            FUN = es_ssgsea, geneset = genes)
        tibble::enframe(ssgsea_vec, name = "sample_id", value = "ssgsea")
    }

    result <- future.apply::future_lapply(
        X = signatures,
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
                    ssgsea = (.data$ssgsea - range(.data$ssgsea)[1]) /
                        (range(.data$ssgsea)[2] - range(.data$ssgsea)[1])
                )
            }
        )
    }
    result <- dplyr::bind_rows(result, .id = "signature_id")
    if (norm == "all") {
        result <- dplyr::mutate(
            result,
            ssgsea = .data$ssgsea /
                (range(.data$ssgsea)[2] - range(.data$ssgsea)[1])
        )
    }

    tidyr::pivot_wider(result,
                       id_cols = "sample_id",
                       names_from = "signature_id",
                       values_from = "ssgsea")

}

.compute_singscore <- function(expr_data, signatures, direction = "none") {
    n_genes <- nrow(expr_data)

    es_singscore <- function(sample_data, geneset) {
        if (direction %in% c("none", "up", "down")) {
            n_in <- length(intersect(geneset, rownames(expr_data)))
            if (direction == "none" | direction == "up") {
                rank_vec <- rank(sample_data)
                if (direction == "none") {
                    rank_center <- ceiling(n_genes / 2)
                    for (i in seq_along(rank_vec)) {
                        rank_vec[[i]] <- abs(rank_vec[[i]] - rank_center)
                    }
                    score <- sum(rank_vec[names(rank_vec) %in% geneset]) / n_in
                    min_score <- (ceiling(n_in / 2) + 1) / 2
                    max_score <- (n_genes - ceiling(n_in / 2) + 1) / 2
                    (score - min_score) / (max_score - min_score)
                }
            } else if (direction == "down") {
                rank_vec <- rank(-sample_data)
            }
            score <- sum(rank_vec[names(rank_vec) %in% geneset]) / n_in
            min_score <- (n_in + 1) / 2
            max_score <- (2 * n_genes - n_in + 1) / 2
            (score - min_score) / (max_score - min_score)
        } else if (direction == "both") {
            rank_up <- rank(sample_data)
            rank_down <- rank(-sample_data)
            min_up <- (n_up + 1) / 2
            max_up <- (2 * n_genes - n_up + 1) / 2
            min_down <- (n_down + 1) / 2
            max_down <- (2 * n_genes - n_down + 1) / 2
            score_up <- sum(rank_up[names(rank_up) %in% geneset]) / n_up
            score_down <- sum(rank_down[names(rank_down) %in% geneset]) / n_down
            score_up <- (score_up - min_up) / (max_up - min_up)
            score_down <- (score_down - min_down) / (max_down - min_down)
            score_up + score_down
        }

        # if (direction == "none") {
        #     rank_vec <- rank(sample_data)
        #     for (i in seq_along(rank_vec)) {
        #         rank_vec[[i]] <- abs(rank_vec[[i]] - rank_center)
        #     }
        #     score <- sum(rank_vec[names(rank_vec) %in% geneset]) / n_in
        #     min_score <- (ceiling(n_in / 2) + 1) / 2
        #     max_score <- (n_genes - ceiling(n_in / 2) + 1) / 2
        #     (score - min_score) / (max_score - min_score)
        # } else if (direction == "up" | direction == "down") {
        #     if (direction == "up") {
        #         rank_vec <- rank(sample_data)
        #     } else {
        #         rank_vec <- rank(-sample_data)
        #     }
        #     score <- sum(rank_vec[names(rank_vec) %in% geneset]) / n_in
        #     min_score <- (n_in + 1) / 2
        #     max_score <- (2 * n_genes - n_in + 1) / 2
        #     (score - min_score) / (max_score - min_score)
        # }


    }

    single_sig_singscore <- function(dataset, genes) {
        singscore_vec <- apply(dataset, MARGIN = 2,
                               FUN = es_singscore, geneset = genes)
        tibble::enframe(singscore_vec, name = "sample_id", value = "singscore")
    }

    result <- future.apply::future_lapply(
        X = signatures,
        FUN = function(genes) {
            single_sig_singscore(dataset = expr_data, genes = genes)
        }
    )
    result <- dplyr::bind_rows(result, .id = "signature_id")

    tidyr::pivot_wider(result,
                       id_cols = "sample_id",
                       names_from = "signature_id",
                       values_from = "singscore")

}
