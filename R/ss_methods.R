#' Single sample scoring methods
#'
#' @description
#' These are internal functions to compute single sample scores from a list of
#'   gene signatures in three different ways:
#'
#'   - **combined z-score** (*Lee et al., 2008*);
#'   - **single sample GSEA** (*Barbie et al., 2009*);
#'   - **singscore** (*Foroutan et al., 2018*).
#'
#'   `compute_ssgsea()` is called by [hack_estimate()] whereas all the three
#'   methods are called by [hack_sig()].
#' @section Algorithm:
#' This section gives a brief explanation of how single sample scores are obtained
#'   from different methods.
#'
#'   ## Combined z-score
#'   Gene expression values are centered by their mean value and scaled by their
#'   standard deviation across samples for each gene (z-scores). Then, for each
#'   sample and signature, corresponding z-scores are added up and divided by the
#'   square root of the signature size (i.e. the number of genes composing a signature).
#'
#'   The combined z-score method is also implemented in the R package `GSVA`
#'   (*Hänzelmann et al., 2013*).
#'
#'   ## Single sample GSEA
#'   For each sample, genes are ranked by expression value in increasing order and
#'   rank normalization may follow (see argument `rank_norm`). Then, two probability-like
#'   vectors are computed for each sample and signature:
#'
#'   - \eqn{P_{in}}, the cumulative sum of weighted ranks divided by their total
#'   sum for genes in the signature;
#'   - \eqn{P_{out}}, the cumulative sum of ones (indicating genes not in the signature)
#'   divided by the number of genes not in the signature.
#'
#'   The single sample GSEA score is obtained by adding up the elements of the
#'   vector difference \eqn{P_{in} - P_{out}}.
#'   Finally, single sample scores could be normalized either across samples or across
#'   gene signatures and samples.
#'
#'   The single sample GSEA method is also implemented in the R package `GSVA`
#'   (*Hänzelmann et al., 2013*).
#'
#'   ## Singscore
#'   For signatures whose genes are supposed to be up- or down-regulated, genes
#'   are ranked by expression value in increasing or decreasing order, respectively.
#'   For signatures whose direction is unknown, genes are ranked by absolute expression
#'   in increasing order and are median-centered.
#'   Enrichment scores are then computed for each sample and signature by averaging
#'   gene ranks for genes in the signature.
#'   Finally, normalized scores are obtained by subtracting the theoretical minimum
#'   mean rank from the score and dividing by the difference between the theoretical
#'   maximum and minimum mean ranks.
#'
#'   The `hacksig` implementation of this method works only with unidirectional (i.e.
#'   all genes up- or down-regulated) and undirected gene signatures.
#'   If you want to get single sample scores for bidirectional gene signatures (i.e.
#'   signatures composed of both up- and down-regulated genes), please use the R
#'   package `singscore` (*Foroutan et al., 2018*).
#' @param signatures A named list of gene signatures.
#' @param direction A character string specifying the **singscore** computation
#'   method depending on the direction of the signatures. Can be on of:
#'
#'   * `"none"` (default), undirected signatures, that is you do not know whether
#'   the genes are up- or down-regulated;
#'   * `"up"`, all genes in the signature are supposed to be up-regulated;
#'   * `"down"`, all genes in the signature are supposed to be down-regulated;
#' @param sample_norm A character string specifying the type of normalization
#'   affecting the **single sample GSEA** scores. Can be one of:
#'
#'   * `"raw"` (default), obtain raw scores;
#'   * `"separate"`, normalize raw scores in \eqn{[0, 1]} across samples for
#'   each signature separately.
#'   * `"all"`, normalize raw scores both across samples and signatures.
#' @param rank_norm A character string specifying how gene expression ranks should
#'   be normalized in the **single sample GSEA** procedure. Valid choices are:
#'
#'   * `"none"` (default), no rank normalization;
#'   * `"rank"`, ranks are multiplied by `10000 / nrow(expr_data)`;
#'   * `"logrank"`, normalized ranks are logged.
#' @param alpha A numeric scalar. Exponent in the running sum of the **single sample GSEA**
#'   score calculation which weighs the gene ranks. Defaults to \eqn{\alpha = 0.25}.
#' @inheritParams hack_estimate
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
#' 108–112. \doi{10.1038/nature08460}.
#'
#' Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., &
#' Davis, M. J. (2018). Single sample scoring of molecular phenotypes.
#' *BMC bioinformatics*, 19(1), 404. \doi{10.1186/s12859-018-2435-4}.
#'
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation
#' analysis for microarray and RNA-seq data. *BMC bioinformatics*, 14, 7.
#' \doi{10.1186/1471-2105-14-7}.
#'
#' Lee, E., Chuang, H. Y., Kim, J. W., Ideker, T., & Lee, D. (2008). Inferring
#' pathway activity toward precise disease classification.
#' *PLoS computational biology*, 4(11), e1000217.
#' \doi{10.1371/journal.pcbi.1000217}.
#'
#' @seealso [hack_sig()], [hack_estimate()], `GSVA::gsva()`, `singscore::multiScore()`
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
        rank_vec <- rank(sample_data)
        if (rank_norm %in% c("rank", "logrank")) {
            rank_vec <- 10000 / n_genes * rank_vec
            if (rank_norm == "logrank") {
                rank_vec <- log(rank_vec + exp(1))
            }
        }
        rank_vec <- rank_vec[order(rank_vec, decreasing = TRUE)]
        prob_in <- rep.int(0, length(rank_vec))
        prob_out <- prob_in
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
            min_score <- (n_in + 1) / 2
            max_score <- (2 * n_genes - n_in + 1) / 2
            if (direction == "none" | direction == "up") {
                rank_vec <- rank(sample_data)
                if (direction == "none") {
                    rank_center <- ceiling(n_genes / 2)
                    for (i in seq_along(rank_vec)) {
                        rank_vec[[i]] <- abs(rank_vec[[i]] - rank_center)
                    }
                    min_score <- (ceiling(n_in / 2) + 1) / 2
                    max_score <- (n_genes - ceiling(n_in / 2) + 1) / 2
                }
            } else if (direction == "down") {
                rank_vec <- rank(-sample_data)
            }
            score <- sum(rank_vec[names(rank_vec) %in% geneset]) / n_in
            (score - min_score) / (max_score - min_score)
        } else stop("Valid choices for 'direction' are 'none', 'up', 'down'",
                    call. = FALSE)
        # else if (direction == "both") {
        #     rank_up <- rank(sample_data)
        #     rank_down <- rank(-sample_data)
        #     min_up <- (n_up + 1) / 2
        #     max_up <- (2 * n_genes - n_up + 1) / 2
        #     min_down <- (n_down + 1) / 2
        #     max_down <- (2 * n_genes - n_down + 1) / 2
        #     score_up <- sum(rank_up[names(rank_up) %in% geneset]) / n_up
        #     score_down <- sum(rank_down[names(rank_down) %in% geneset]) / n_down
        #     score_up <- (score_up - min_up) / (max_up - min_up)
        #     score_down <- (score_down - min_down) / (max_down - min_down)
        #     score_up + score_down
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
