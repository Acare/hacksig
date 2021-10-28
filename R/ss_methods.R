#' Single sample methods
#'
#' @description
#' These are internal functions to compute single sample scores in three different
#' ways: **z-score** (*Lee et al., 2008*), **ssGSEA** (*Barbie et al., 2009*)
#' or **singscore** (*Foroutan et al., 2018*).
#' `compute_ssgsea` is called by [hack_estimate()] whereas all the three
#' methods are called by \code{\link{hack_sig}}.
#'
#' @section Algorithm:
#' ## Z-score:
#' scale and sum z-scores
#' ## Single sample GSEA:
#' ranks and running sum
#' ## Singscore:
#' ranks and mean
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'   row names and samples as columns.
#' @param signatures It can be a list of signatures or a character string indicating
#'   a keyword for a group of signatures (e.g. "immune" or "ifng" for immune and
#'   interferon gamma signatures respectively). The default (`all`) will cause the
#'   function to compute single sample scores for all the signatures implemented
#'   in `hacksig`.
#' @param direction A character specifying the **singscore** computation depending on the
#'   direction of the signatures. Can be on of:
#'
#'   * `none`, undirected signatures, that is you don't know whether the genes are
#'     up- or down-regulated (default);
#'   * `up`, all genes in the signature are supposed to be up-regulated;
#'   * `down`, all genes in the signature are supposed to be down-regulated;
#'   * `both`, a signature is composed of up- and down sub-signatures.
#'     You must supply it as a nested list.
#' @param sample_norm A character string specifying the type of normalization affecting
#'   the **single sample GSEA** scores. Can be one of:
#'
#'   * `raw`, obtain raw scores (default);
#'   * `separate`, normalize raw scores in \eqn{[0, 1]} across samples for each signature separately.
#'   * `all`, normalize raw scores both across samples and signatures.
#' @param rank_norm A character string specifying how gene expression ranks should
#'   be normalized. Valid choices are:
#'
#'   * `none`, no rank normalization (default);
#'   * `rank`, ranks are multiplied by `10000 / nrow(expr_data)`;
#'   * `logrank`,
#' @param alpha A numeric scalar. Exponent in the running sum of the **single sample GSEA**
#'   score calculation which weights the gene ranks. Defaults to \eqn{\alpha = 0.25}.
#'
#' @return A tibble with one row for each sample in `expr_data`, a column `sample_id`
#'   indicating sample identifiers and one column for each input signature giving
#'   single sample scores.
#'
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
#'
#' @seealso [hack_sig()], [hack_estimate()]
#'
#' @keywords internal
#' @name ss_methods
NULL

#' @rdname ss_methods
compute_zscore <- function(expr_data, signatures) {
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
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

#' @rdname ss_methods
compute_ssgsea <- function(expr_data, signatures, sample_norm = "raw",
                           rank_norm = "none", alpha = 0.25) {
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
    n_genes <- nrow(expr_data)

    es_ssgsea <- function(sample_data, geneset) {
        n_in <- length(intersect(geneset, rownames(expr_data)))
        n_out <- n_genes - n_in
        # check rank options as in GSVA::gsva()
        # rank_vec <- rank(sample_data, na.last = "keep", ties.method = "average")
        rank_vec <- rank(sample_data)
        if (rank_norm %in% c("rank", "logrank")) {
            rank_vec <- 10000 / n_genes * rank_vec
            if (rank_norm == "logrank") {
                rank_vec <- log(rank_vec + exp(1))
            }
        }
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
    if (sample_norm == "separate") {
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
    if (sample_norm == "all") {
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

#' @rdname ss_methods
compute_singscore <- function(expr_data, signatures, direction = "none") {
    if (is.matrix(expr_data) == TRUE) {
        expr_data <- as.data.frame(expr_data)
    }
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
