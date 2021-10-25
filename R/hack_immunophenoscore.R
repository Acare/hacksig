#' Immunophenoscore
#'
#' Obtain various immune biomarkers scores, which combined together give the
#'  immunophenoscore (*Charoentong et al., 2017*).
#'
#' @inheritParams hack_estimate
#' @return A data frame with
#'
#' @references
#' Charoentong, P., Finotello, F., Angelova, M., Mayer, C.,
#' Efremova, M., Rieder, D., Hackl, H., & Trajanoski, Z. (2017). Pan-cancer
#' Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and
#' Predictors of Response to Checkpoint Blockade. *Cell reports*, 18(1), 248–262.
#' [doi: 10.1016/j.celrep.2016.12.019](https://doi.org/10.1016/j.celrep.2016.12.019).
#'
#' @importFrom rlang .data
#' @export
hack_immunophenoscore <- function(expr_data) {
    signatures_data <- hacksig::signatures_data
    biom_classes <- tibble::tibble(
        gene_type = c("B2M", "TAP1", "TAP2",
                      paste0("HLA-", c(LETTERS[1:3], "DPA1", "DPB1", "E", "F")),
                      "PD-1", "CTLA-4", "LAG3", "TIGIT", "TIM3", "PD-L1", "PD-L2",
                      "CD27", "ICOS", "IDO1",
                      "Act CD4", "Act CD8", "Tem CD4", "Tem CD8", "MDSC", "Treg"),
        gene_class = c(rep(c("MHC", "CP"), each = 10),
                       rep_len("EC", 4), "SC", "SC")
    )
    ips_genes <- signatures_data[signatures_data$signature_id == "immunophenoscore",
                                 c("gene_type", "gene_symbol", "gene_weight")]
    ips_genes <- merge(ips_genes, biom_classes, by = "gene_type")
    scaled_data <- scale(expr_data, center = TRUE, scale = TRUE)
    ips_data <- merge(ips_genes,
                      tibble::as_tibble(scaled_data, rownames = "gene_symbol"),
                      by = "gene_symbol")
    # browser()
    ips_data <- tidyr::pivot_longer(
        ips_data,
        cols = -c("gene_symbol", "gene_type", "gene_weight", "gene_class"),
        names_to = "sample_id",
        values_to = "norm_expr"
        )
    ips_data <- dplyr::mutate(dplyr::group_by(ips_data, .data$sample_id, .data$gene_type),
                              type_score = mean(.data$norm_expr),
                              type_weight = mean(.data$gene_weight))
    keep_cols <- c("sample_id", "gene_type", "gene_class", "type_score", "type_weight")
    result_type <- dplyr::distinct(dplyr::ungroup(ips_data[, keep_cols]))
    result_type$weighted_type_score <- result_type$type_weight * result_type$type_score
    result_type$type_weight <- NULL
    result_class <- dplyr::mutate(dplyr::group_by(result_type, .data$sample_id, .data$gene_class),
                                  class_score = mean(.data$weighted_type_score))
    keep_cols <- c("sample_id", "gene_class", "class_score")
    result_class <- dplyr::distinct(dplyr::ungroup(result_class[, keep_cols]))
    result_class <- dplyr::mutate(
        dplyr::group_by(result_class, .data$sample_id),
        raw_score = sum(.data$class_score),
        ips_score = dplyr::case_when(
            .data$raw_score <= 0 ~ 0,
            .data$raw_score >= 3 ~ 10,
            dplyr::between(.data$raw_score, 0, 3) ~ round(.data$raw_score * 10 / 3, digits = 0)
            )
        )
    dplyr::left_join(result_type, result_class, by = c("sample_id", "gene_class"))
}
