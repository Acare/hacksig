#' CINSARC
#'
#' Given a gene expression data frame and a vector indicating the distant metastasis
#' status of a sample, `hack_cinsarc()` classifies samples in one of two risk
#' classes, C1 or C2, using the 67-gene signature CINSARC and the LOOCV method
#' as implemented in *Chibon et al. (2010), Lesluyes & Chibon (2020)*.
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param dm_status A numeric vector specifying whether a sample has either (1)
#'  or not (0) developed distant metastasis.
#'
#' @return A data frame with
#' @references
#' Chibon, F., Lagarde, P., Salas, S., Pérot, G., Brouste, V., Tirode, F.,
#' Lucchesi, C., de Reynies, A., Kauffmann, A., Bui, B., Terrier, P.,
#' Bonvalot, S., Le Cesne, A., Vince-Ranchère, D., Blay, J. Y., Collin, F.,
#' Guillou, L., Leroux, A., Coindre, J. M., & Aurias, A. (2010). Validated
#' prediction of clinical outcome in sarcomas and multiple types of cancer on
#' the basis of a gene expression signature related to genome complexity.
#' *Nature medicine*, 16(7), 781–787.
#' [doi: 10.1038/nm.2174](https://doi.org/10.1038/nm.2174).
#'
#' Lesluyes, T., & Chibon, F. (2020). A Global and Integrated
#' Analysis of CINSARC-Associated Genetic Defects. *Cancer research*, 80(23),
#' 5282–5290. [doi: 10.1158/0008-5472.CAN-20-0512](https://doi.org/10.1158/0008-5472.CAN-20-0512).
#'
#' @examples
#' hack_cinsarc
#'
#' @export
hack_cinsarc <- function(expr_data, dm_status) {
    signatures_data <- hacksig::signatures_data
    event_df <- tibble::tibble(sample_id = colnames(expr_data),
                               event = as.factor(ifelse(dm_status == 0, "X0", "X1")))
    cinsarc_genes <- signatures_data[signatures_data$signature_id == "cinsarc",
                                     "gene_symbol", drop = TRUE]
    filt_data <- expr_data[rownames(expr_data) %in% cinsarc_genes, ]
    centered_data <- scale(t(filt_data), center = TRUE, scale = FALSE)
    result <- vector("list", ncol(expr_data))

    for (i in seq_along(event_df$sample_id)) {
        loo_scale <- scale(t(filt_data[, -i]), center = TRUE, scale = FALSE)
        loo_scale <- tibble::as_tibble(loo_scale, rownames = "sample_id")
        loo_scale <- merge(loo_scale, event_df, by = "sample_id")
        browser()
        loo_x0 <- loo_scale[loo_scale$event == "X0",
                            unlist(lapply(loo_scale, is.numeric), use.names = FALSE)]
        loo_x1 <- loo_scale[loo_scale$event == "X1",
                            unlist(lapply(loo_scale, is.numeric), use.names = FALSE)]
        loo_centroids <- data.frame(x0 = colMeans(loo_x0),
                                    x1 = colMeans(loo_x1))
        loo_cor <- 1 - stats::cor(t(centered_data), loo_centroids,
                                  method = "spearman")
        loo_cor <- tibble::as_tibble(loo_cor, rownames = "sample_id")
        loo_cor$cinsarc_class <- ifelse(loo_cor$x0 < loo_cor$x1, "C1", "C2")
        result[[i]] <- loo_cor[i, ]
    }

    dplyr::bind_rows(result)

}
