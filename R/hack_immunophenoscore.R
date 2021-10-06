#' Classify a cohort of samples based on the Immunophenoscore
#'
#' WRITE SOMETHING
#'
#' @param expr_data A gene expression matrix (or data frame) with gene symbols as
#'  row names and samples as columns.
#'
#' @return A data frame with
#' @export
#'
#' @examples
#' m <- matrix()
hack_immunophenoscore <- function(expr_data) {
    IPSG <- read_tsv(file = path_to_IPS_genes)
    expr_mat <- expression_matrix %>% as.data.frame() %>% rownames_to_column("GENE")

    # compute sample-wise mean and sd for all genes in expression_matrix
    mean_sd <- expression_matrix %>%
        t() %>%
        as_tibble(rownames = "sample") %>%
        rowwise() %>%
        mutate(mean = mean(c_across(-sample)), sd = sd(c_across(-sample)), .keep = "unused")

    # check if all genes are present and eventually print out a message
    # remove missing immune genes
    MISSING_GENES <- setdiff(IPSG$GENE, expr_mat$GENE)

    if (length(MISSING_GENES) > 0) {
        IPSG <- IPSG %>% filter(!GENE %in% MISSING_GENES)
        cat("differently named or missing genes: ", sort(MISSING_GENES), "\n")
    }

    # filter immune genes and compute z-scores
    immune_mat <- IPSG %>%
        left_join(expr_mat, by = "GENE") %>%
        pivot_longer(-c(GENE, NAME, CLASS, WEIGHT), names_to = "sample", values_to = "expr") %>%
        left_join(mean_sd, by = "sample") %>%
        mutate(z_score = (expr - mean) / sd)

    # create summary data frame
    IPS_df <- immune_mat %>%
        group_by(sample, NAME) %>%
        mutate(score_NAME = mean(z_score), weight_NAME = mean(WEIGHT)) %>%
        ungroup() %>%
        select(-GENE, -expr, -mean, -sd, -z_score) %>%
        distinct() %>%
        # filter(sample == "TCGA-04-1348") %>%
        mutate(score_NAME_w = score_NAME * weight_NAME) %>%
        group_by(sample, CLASS) %>%
        mutate(score_CLASS = mean(score_NAME_w)) %>%
        group_by(sample) %>%
        mutate(sum_CLASS = sum(unique(score_CLASS)),
               IPS = case_when(sum_CLASS <= 0 ~ 0,
                               sum_CLASS >= 3 ~ 10,
                               sum_CLASS > 0 | sum_CLASS < 3 ~ round(sum_CLASS * 10 / 3, digits = 0))) %>%
        select(sample, NAME, CLASS, starts_with("score"), IPS) %>%
        ungroup()

    IPS_df
}
