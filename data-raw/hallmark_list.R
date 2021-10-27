library(tidyverse)

hallmark_list <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct() %>%
    nest(genes = c(gene_symbol)) %>%
    mutate(genes = map(genes, compose(as_vector, unname))) %>%
    deframe()

usethis::use_data(hallmark_list, overwrite = TRUE)
