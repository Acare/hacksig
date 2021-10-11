library(dplyr)
library(org.Hs.eg.db)

signatures_data <- readr::read_csv("~/Desktop/hacksig_ideas/temp_signatures.csv",
                                   col_types = "ccccnc")
# mutate(gene_symbol = case_when(
#     gene_symbol == "WISP1" ~ "CCN4",
#     gene_symbol == "GPR124" ~ "ADGRA2",
#     gene_symbol == "LPPR4" ~ "PLPPR4",
#     gene_symbol == "TXNDC3" ~ "NME8",
#     gene_symbol == "ODZ4" ~ "TENM4",
#     gene_symbol == "FYB" ~ "FYB1",
#     TRUE ~ gene_symbol
# ))

# add entrez ids
query <- select(org.Hs.eg.db,
                keys = unique(tab_sig$gene_symbol),
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL")

signatures_data <- signatures_data %>%
    left_join(query, by = c("gene_symbol" = "SYMBOL")) %>%
    mutate(entrez_gene_id = ENTREZID, .keep = "unused")

signatures_data %>% readr::write_csv("~/Desktop/hacksig_ideas/temp_signatures.csv")

usethis::use_data(signatures_data, overwrite = TRUE)