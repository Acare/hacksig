#' Hack the CINSARC classification
#'
#' @description
#' Given a gene expression matrix and a 0-1 vector indicating the distant metastasis
#'   status of samples, `hack_cinsarc()` classifies samples into one of two risk
#'   classes, C1 or C2, using the CINSARC signature as implemented in
#'   *Chibon et al., 2010*.
#' @details
#' CINSARC (*Complexity INdex in SARComas*) is a prognostic 67-gene signature
#'   related to mitosis and control of chromosome integrity.
#'   It was developed to improve metastatic outcome prediction in soft tissue
#'   sarcomas over the FNCLCC (*Fédération Francaise des Centres de Lutte Contre
#'   le Cancer*) grading system.
#' @section Algorithm:
#' The CINSARC method implemented in `hacksig` makes use of leave-one-out cross
#'   validation (LOOCV) to classify samples into C1/C2 risk groups (see *Lesluyes & Chibon, 2020*).
#'   First, gene expression values are centered by their mean across samples.
#'   Then, for each iteration of the LOOCV, mean normalized gene values are computed
#'   by metastasis group (i.e. compute the metastatic centroids). Then, one minus the
#'   Spearman's correlation between centered samples and metastatic centroids are computed.
#'   Finally, if a sample is more correlated to the non-metastatic centroid, then
#'   it is assigned to the C1 class (low risk). Conversely, if a sample is more
#'   correlated to the metastatic centroid, then it is assigned to the C2 class (high risk).
#' @param dm_status A numeric vector specifying whether a sample has either (1)
#'   or not (0) developed distant metastasis.
#' @inheritParams hack_estimate
#' @return A tibble with one row for each sample in `expr_data` and two columns:
#'   `sample_id` and `cinsarc_class`.
#' @source [codeocean.com/capsule/4933686/tree/v4](https://codeocean.com/capsule/4933686/tree/v4)
#'
#' @references
#' Chibon, F., Lagarde, P., Salas, S., Pérot, G., Brouste, V., Tirode, F.,
#' Lucchesi, C., de Reynies, A., Kauffmann, A., Bui, B., Terrier, P.,
#' Bonvalot, S., Le Cesne, A., Vince-Ranchère, D., Blay, J. Y., Collin, F.,
#' Guillou, L., Leroux, A., Coindre, J. M., & Aurias, A. (2010). Validated
#' prediction of clinical outcome in sarcomas and multiple types of cancer on
#' the basis of a gene expression signature related to genome complexity.
#' *Nature medicine*, 16(7), 781–787. \doi{10.1038/nm.2174}.
#'
#' Lesluyes, T., & Chibon, F. (2020). A Global and Integrated
#' Analysis of CINSARC-Associated Genetic Defects. *Cancer research*, 80(23),
#' 5282–5290. \doi{10.1158/0008-5472.CAN-20-0512}.
#'
#' @examples
#' # generate random distant metastasis outcome
#' set.seed(123)
#' test_dm_status <- sample(c(0, 1), size = ncol(test_expr), replace = TRUE)
#'
#' hack_cinsarc(test_expr, test_dm_status)
#' @export
hack_cinsarc <- function(expr_data, dm_status) {
    sig_data <- signatures_data
    event_df <- tibble::tibble(sample_id = colnames(expr_data),
                               event = as.factor(ifelse(dm_status == 0, "X0", "X1")))
    cinsarc_genes <- sig_data[sig_data$signature_id == "cinsarc",
                              "gene_symbol", drop = TRUE]
    filt_data <- expr_data[rownames(expr_data) %in% cinsarc_genes, ]
    centered_data <- scale(t(filt_data), center = TRUE, scale = FALSE)
    result <- vector("list", ncol(expr_data))
    for (i in seq_along(event_df$sample_id)) {
        loo_scale <- scale(t(filt_data[, -i]), center = TRUE, scale = FALSE)
        loo_scale <- tibble::as_tibble(loo_scale, rownames = "sample_id")
        loo_scale <- merge(loo_scale, event_df, by = "sample_id")
        loo_x0 <- loo_scale[loo_scale$event == "X0",
                            unlist(lapply(loo_scale, is.numeric), use.names = FALSE)]
        loo_x1 <- loo_scale[loo_scale$event == "X1",
                            unlist(lapply(loo_scale, is.numeric), use.names = FALSE)]
        loo_centroids <- data.frame(x0 = colMeans(loo_x0), x1 = colMeans(loo_x1))
        loo_cor <- 1 - stats::cor(t(centered_data), loo_centroids,
                                  method = "spearman")
        loo_cor <- tibble::as_tibble(loo_cor, rownames = "sample_id")
        loo_cor$cinsarc_class <- ifelse(loo_cor$x0 < loo_cor$x1, "C1", "C2")
        result[[i]] <- loo_cor[i, c("sample_id", "cinsarc_class")]
    }
    dplyr::bind_rows(result)
}
