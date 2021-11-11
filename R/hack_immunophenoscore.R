#' Hack the Immunophenoscore
#'
#' @description
#' Obtain various immune biomarkers scores, which combined together give the
#'  immunophenoscore (*Charoentong et al., 2017*).
#'
#' @details
#' Immunophenoscore (IPS), is a score that was generated, through the use of machine learning, with The Cancer Genome Atlas data
#' from 20 solid cancers. The Cancer Immunome Atlas was generated (https://tcia.at/) by the depiction of intratumoral immune landscapes.
#' The score identifies the tumor immunogenicity, and it was able to predict the response to immunocheckpoint inhibitors
#' (anti-cytotoxic T lymphocyte antigen-4 (CTLA-4) and anti-programmed cell death protein 1 (anti-PD-1) antibodies) in two validation cohorts.
#'
#' @param feature A string indicating whether you want a more granular output (`type`)
#'   or a more aggregated one (`class`). Each of the two possible choices will give
#'   the raw and discrete immunophenoscores.
#' @inheritParams hack_estimate
#' @return A data frame with
#'
#' @source [github.com/icbi-lab/Immunophenogram](https://github.com/icbi-lab/Immunophenogram)
#'
#' @references
#' Charoentong, P., Finotello, F., Angelova, M., Mayer, C.,
#' Efremova, M., Rieder, D., Hackl, H., & Trajanoski, Z. (2017). Pan-cancer
#' Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and
#' Predictors of Response to Checkpoint Blockade. *Cell reports*, 18(1), 248–262.
#' [doi: 10.1016/j.celrep.2016.12.019](https://doi.org/10.1016/j.celrep.2016.12.019).
#'
#' @importFrom rlang .data
#' @seealso [hack_sig()] to compute Immunophenoscore biomarkers in different
#'   ways (use `signatures = "ips"`).
#'
#'   [check_sig()] to check if all/most of the Immunophenoscore biomarkers are
#'   present in your expression matrix (use `signatures = "ips"`).
#' @export
hack_immunophenoscore <- function(expr_data, feature = "all") {
    sig_data <- hacksig::signatures_data
    biom_classes <- tibble::tibble(
        gene_type = c("b2m", "tap1", "tap2",
                      paste0("hla_", c(letters[1:3], "dpa1", "dpb1", "e", "f")),
                      "pd_1", "ctla_4", "lag3", "tigit", "tim3", "pd_l1", "pd_l2",
                      "cd27", "icos", "ido1",
                      "act_cd4", "act_cd8", "tem_cd4", "tem_cd8", "mdsc", "treg"),
        gene_class = c(rep(c("mhc", "cp"), each = 10),
                       rep_len("ec", 4), "sc", "sc")
    )
    ips_genes <- sig_data[grep("ips", sig_data$signature_id),
                          c("signature_id", "gene_symbol", "gene_weight")]
    ips_genes$signature_id <- gsub("ips_", "", ips_genes$signature_id)
    names(ips_genes)[names(ips_genes) == "signature_id"] <- "gene_type"
    ips_genes <- dplyr::left_join(ips_genes, biom_classes, by = "gene_type")
    scaled_data <- scale(expr_data, center = TRUE, scale = TRUE)
    ips_data <- merge(ips_genes,
                      tibble::as_tibble(scaled_data, rownames = "gene_symbol"),
                      by = "gene_symbol")
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
        raw = sum(.data$class_score),
        ips = dplyr::case_when(
            .data$raw <= 0 ~ 0,
            .data$raw >= 3 ~ 10,
            .data$raw > 0 | .data$raw < 3 ~ round(.data$raw * 10 / 3, digits = 0)
            )
        )
    result_class <- dplyr::ungroup(
        tidyr::pivot_wider(
            result_class,
            c("sample_id", "raw", "ips"),
            names_from = "gene_class",
            values_from = "class_score"
        )
    )
    names(result_class)[-1] <- paste0(names(result_class)[-1], "_score")
    names(result_class) <- tolower(names(result_class))
    result_type <- tidyr::pivot_wider(
        result_type,
        "sample_id",
        names_from = "gene_type",
        values_from = "weighted_type_score"
    )
    names(result_type) <- tolower(gsub(" |-", "_", names(result_type)))
    names(result_type)[-1] <- paste0(names(result_type)[-1], "_score")
    result <- dplyr::left_join(result_class, result_type, by = "sample_id")
    if (feature == "class") {
        result_class
    } else if (feature == "type") {
        result[, -c("ec_score", "mhc_score", "sc_score", "cp_score")]
    } else {
        result
    }
}
