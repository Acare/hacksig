#' Gene-level information for each signature implemented
#'
#' A data set used by the `hacksig` functions to retrieve gene-level information
#' about XXX signatures implemented in the package.
#'
#' @format A tibble with six columns:
#'
#'   * `signature_id`: a unique ID given to the signature;
#'   * `signature_keyword`: a list of keywords which can be used to filter signatures;
#'   * `gene_symbol`: HUGO Gene Symbol ID;
#'   * `gene_entrez_id`: Entrez gene ID;
#'   * `gene_weight`: gene weights used in the signature model;
#'   * `signature_method`: default method implemented in the original publication.
#' @details
#' Identifiers in `signature_id` correspond to the following DOIs:
#'
#'   * `"rooney2015_cytoact"`, [10.1016/j.cell.2014.12.033](https://doi.org/10.1016/j.cell.2014.12.033)
#'   * `"ayers2017_immexp"`, [10.1172/JCI91190](https://doi.org/10.1172/JCI91190)
#'   * `"she2020_irgs"`, [10.1186/s12935-020-1104-7](https://doi.org/10.1186/s12935-020-1104-7)
#'   * `"wu2020_metabolic"`, [10.1155/2020/6716908](https://doi.org/10.1155/2020/6716908)
#'   * `"fang2020_irgs"`, [10.18632/aging.202842](https://doi.org/10.18632/aging.202842)
#'
#' @seealso [hack_sig()], [check_sig()]
#' @docType data
#' @keywords datasets
"signatures_data"


#' A toy gene expression matrix
#'
#' A gene expression matrix simulating expression profiles of 20 samples. It should
#' be used for testing purpose.
#'
#' @format A random normal data matrix with 20000 genes as rows and 20 samples as columns.
#' @examples
#' hack_sig(test_expr)
#' @docType data
#' @keywords datasets
"test_expr"

#' Hallmark gene set list
#'
#' A gene expression matrix ...
#'
#' @format A named list storing the 50 Hallmark gene sets with gene symbols as
#'   gene identifiers.
#' @examples
#' hack_sig(test_expr, hallmark_list, method = "zscore")
#' @docType data
#' @keywords datasets
#' @source \url{http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=H}
"hallmark_list"
