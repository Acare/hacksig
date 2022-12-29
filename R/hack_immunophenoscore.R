#' Hack the Immunophenoscore
#'
#' @description
#' Obtain various immune biomarkers scores, which combined together give the
#'   immunophenoscore (*Charoentong et al., 2017*).
#' @details
#' The immunophenoscore is conceived as a quantification of tumor immunogenicity.
#'   It is obtained by aggregating multiple immune biomarkers scores, which are
#'   grouped into four major classes:
#'
#'   * _MHC molecules_ (__MHC__), expression of MHC class I, class II, and non-classical molecules;
#'   * _Immunomodulators_ (__CP__), expression of certain co-inhibitory and co-stimulatory molecules;
#'   * _Effector cells_ (__EC__), infiltration of activated CD8+/CD4+ T cells and Tem (effector memory) CD8+/CD4+ cells;
#'   * _Suppressor cells_ (__SC__), infiltration of immunosuppressive cells (Tregs and MDSCs).
#'
#'   The table below shows in detail the 26 immune biomarkers and cell types grouped
#'   by class together with the number of genes which represent them:
#'
#'   | __Class__ \|| __Biomarker/cell type__ \|| __No. genes__ |
#'   | ----- |:-------------------:|:---------:|
#'   | MHC   | B2M                 | 1         |
#'   | MHC   | HLA-A               | 1         |
#'   | MHC   | HLA-B               | 1         |
#'   | MHC   | HLA-C               | 1         |
#'   | MHC   | HLA-DPA1            | 1         |
#'   | MHC   | HLA-DPB1            | 1         |
#'   | MHC   | HLA-E               | 1         |
#'   | MHC   | HLA-F               | 1         |
#'   | MHC   | TAP1                | 1         |
#'   | MHC   | TAP2                | 1         |
#'   | CP    | CD27                | 1         |
#'   | CP    | CTLA-4              | 1         |
#'   | CP    | ICOS                | 1         |
#'   | CP    | IDO1                | 1         |
#'   | CP    | LAG3                | 1         |
#'   | CP    | PD1                 | 1         |
#'   | CP    | PD-L1               | 1         |
#'   | CP    | PD-L2               | 1         |
#'   | CP    | TIGIT               | 1         |
#'   | CP    | TIM3                | 1         |
#'   | EC    | Act CD4             | 24        |
#'   | EC    | Act CD8             | 26        |
#'   | EC    | Tem CD4             | 27        |
#'   | EC    | Tem CD8             | 25        |
#'   | SC    | MDSC                | 20        |
#'   | SC    | Treg                | 20        |
#' @section Algorithm:
#' Samplewise gene expression z-scores are obtained for each of 26 immune cell
#'   types and biomarkers. Then, weighted averaged z-scores are computed for each
#'   class and the raw immunophenoscore (\eqn{IPS_{raw}}) results as the sum of the
#'   four class scores. Finally, the immunophenoscore (\eqn{IPS}) is given as an
#'   integer value between 0 and 10 in the following way:
#'
#'   - \eqn{IPS = 0}, if \eqn{IPS_{raw} \le 0};
#'   - \eqn{IPS = [10 * (IPS_{raw} / 3)]}, if \eqn{0 < IPS_{raw} < 3};
#'   - \eqn{IPS = 10}, if \eqn{IPS_{raw} \ge 3}.
#' @param extract A string controlling which type of biomarker scores you want
#'   to obtain. Possible choices are:
#'
#'   - `"ips"` (default), only raw and discrete IPS scores;
#'   - `"class"`, IPS scores together with the four summary class scores;
#'   - `"all"`, all possible biomarker scores.
#' @inheritParams hack_estimate
#' @return A tibble with one row for each sample in `expr_data`, a column `sample_id`
#'   indicating sample identifiers and a number of additional columns depending
#'   on the choice of `extract`.
#' @examples
#' hack_immunophenoscore(test_expr)
#' hack_immunophenoscore(test_expr, extract = "class")
#' @source [github.com/icbi-lab/Immunophenogram](https://github.com/icbi-lab/Immunophenogram)
#' @references
#' Charoentong, P., Finotello, F., Angelova, M., Mayer, C.,
#' Efremova, M., Rieder, D., Hackl, H., & Trajanoski, Z. (2017). Pan-cancer
#' Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and
#' Predictors of Response to Checkpoint Blockade. *Cell reports*, 18(1), 248–262.
#' \doi{10.1016/j.celrep.2016.12.019}.
#' @seealso [hack_sig()] to compute Immunophenoscore biomarkers in different
#'   ways (e.g. use `signatures = "ips"` and `method = "singscore"`).
#'
#'   [check_sig()] to check if all/most of the Immunophenoscore biomarkers are
#'   present in your expression matrix (use `signatures = "ips"`).
#' @importFrom data.table `:=`
#' @export
hack_immunophenoscore <- function(expr_data, extract = "ips") {
    sig_data <- signatures_data
    class_score = gene_weight = norm_expr = weighted_type_score = NULL # due to NSE notes in R CMD check
    ips_genes <- sig_data[grep("ips", sig_data$signature_id),
                          c("signature_keywords", "signature_id", "gene_symbol", "gene_weight")]
    ips_genes$signature_keywords <- gsub("\\|.*", "", ips_genes$signature_keywords)
    ips_genes$signature_keywords <- gsub(" |-", "_", ips_genes$signature_keywords)
    ips_genes$signature_id <- gsub("ips_", "", ips_genes$signature_id)
    names(ips_genes)[names(ips_genes) == "signature_keywords"] <- "gene_type"
    names(ips_genes)[names(ips_genes) == "signature_id"] <- "gene_class"
    scaled_data <- scale(expr_data, center = TRUE, scale = TRUE)
    ips_data <- merge(ips_genes,
                      tibble::as_tibble(scaled_data, rownames = "gene_symbol"),
                      by = "gene_symbol")
    data.table::setDT(ips_data)
    ips_data <- data.table::melt(ips_data,
                                 id.vars = c("gene_symbol", "gene_type", "gene_weight", "gene_class"),
                                 variable.name = "sample_id",
                                 value.name = "norm_expr")
    ips_data[,
             c("type_score", "type_weight") :=
                 list(mean(norm_expr, na.rm = TRUE), mean(gene_weight, na.rm = TRUE)),
             by = c("sample_id", "gene_type")]
    result_type <- unique(
        ips_data[, c("sample_id", "gene_type", "gene_class", "type_score", "type_weight")]
    )
    result_type$weighted_type_score <- result_type$type_weight * result_type$type_score
    result_type[, "type_weight" := NULL]

    result_class <- data.table::copy(result_type)
    result_class[,
                 "class_score" := mean(weighted_type_score, na.rm = TRUE),
                 by = c("sample_id", "gene_class")]
    result_class <- unique(
        result_class[, c("sample_id", "gene_class", "class_score")]
    )
    result_class[
        ,
        "raw" := sum(class_score, na.rm = TRUE),
        by = "sample_id"
    ][
        ,
        "ips" := data.table::fcase(
            raw <= 0, 0,
            raw >= 3, 10,
            raw > 0 | raw < 3, round(raw * 10 / 3, digits = 0)
        ),
        by = "sample_id"
    ]

    result_class <- data.table::dcast(
            result_class,
            sample_id + raw + ips ~ gene_class,
            value.var = "class_score"
    )
    result_type <- data.table::dcast(
        result_type,
        sample_id ~ gene_type,
        value.var = "weighted_type_score"
    )
    names(result_class)[-1] <- paste0(names(result_class)[-1], "_score")
    names(result_class) <- tolower(names(result_class))
    names(result_type) <- tolower(gsub(" |-", "_", names(result_type)))
    names(result_type)[-1] <- paste0(names(result_type)[-1], "_score")

    result <- result_type[result_class, on = "sample_id"]

    if (extract == "class") {
        tibble::as_tibble(result_class)
    } else if (extract == "ips") {
        tibble::as_tibble(result_class[, c("sample_id", "raw_score", "ips_score")])
    } else if (extract == "all") {
        tibble::as_tibble(result)
    } else stop("Must provide a valid string for 'extract'.
                Possible choices are 'ips', 'class', 'all'.", call. = FALSE)
}
