#' Classify a cohort of samples based on the CINSARC Signature
#'
#' Given a gene expression data frame and a vector indicating the distant metastasis
#' status of a sample, `hack_cinsarc()` classifies samples in one of two risk
#' classes, C1 or C2, using the 67-gene signature CINSARC and the LOOCV method
#' as implemented in ...
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#' @param dm_status A numeric vector specifying whether a sample has either (1)
#'  or not (0) developed distant metastasis.
#'
#' @return A data frame with
#' @export
#'
#' @examples
#' m <- matrix()
hack_cinsarc <- function(expr_data, dm_status) {
    event_df <- tibble::tibble(sample_name = colnames(expr_data),
                               event = ifelse(event_vector == 0, "X0", "X1") %>% as.factor())

    tidy_transpose <- function(data, genes) {
        data %>%
            as.data.frame() %>%
            tibble::rownames_to_column("SYMBOL") %>%
            dplyr::filter(SYMBOL %in% genes) %>%
            tidyr::pivot_longer(-SYMBOL, names_to = "sample_name") %>%
            tidyr::pivot_wider(names_from = SYMBOL, values_from = value)
    }

    expr_mat_centered <- expression_matrix %>%
        tidy_transpose(CINSARC) %>%
        dplyr::mutate(across(where(is.numeric), ~ .x - mean(.x))) %>%
        tidyr::pivot_longer(-sample_name, names_to = "SYMBOL") %>%
        tidyr::pivot_wider(names_from = sample_name, values_from = value) %>%
        tibble::column_to_rownames("SYMBOL")

    missing_genes <- setdiff(CINSARC, rownames(expr_mat_centered))

    if (length(missing_genes) > 0) {
        message("The following CINSARC genes are missing:\n", stringr::str_c(missing_genes, collapse = ", "))
    }

    class_pred <- tibble(sample_name = colnames(expr_mat_centered),
                         X0 = NA_real_, X1 = NA_real_, CINSARC_pred = NA_character_)

    for (i in 1:ncol(expr_mat_centered)) {

        one_out_centers <- expression_matrix %>%
            tidy_transpose(CINSARC) %>%
            dplyr::slice(-i) %>%
            dplyr::mutate(across(where(is.numeric), ~ .x - mean(.x))) %>%
            dplyr::left_join(event_df, by = "sample_name") %>%
            dplyr::group_by(event) %>%
            dplyr::summarize(across(where(is.numeric), mean)) %>%
            tidyr::pivot_longer(-event, names_to = "SYMBOL") %>%
            tidyr::pivot_wider(names_from = event, values_from = value) %>%
            tibble::column_to_rownames("SYMBOL")

        pred_vec <- (1 - cor(expr_mat_centered, one_out_centers, method = "spearman")) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("sample_name") %>%
            dplyr::rowwise() %>%
            dplyr::mutate(CINSARC_pred = ifelse(X0 < X1, "C1", "C2")) %>%
            dplyr::ungroup() %>%
            dplyr::slice(i)

        class_pred <- class_pred %>% dplyr::rows_update(pred_vec, by = "sample_name")

    }

    class_pred
}
