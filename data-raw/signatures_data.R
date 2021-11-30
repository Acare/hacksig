library(dplyr)
library(org.Hs.eg.db) # version 3.12.0
# library(hgu133plus2.db) # version 3.2.3

# Add some signatures to hacksig_signatures.csv. Then:
signatures_data <- readr::read_csv("data-raw/hacksig_signatures.csv",
                                   col_types = "ccccnccc")

genes_to_update <- signatures_data %>%
    filter(is.na(gene_entrez_id)) %>%
    pull(gene_symbol) %>%
    unique()

# Add entrez ids from gene symbols
query <- select(org.Hs.eg.db,
                keys = genes_to_update,
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL") %>%
    rename(gene_symbol = SYMBOL, gene_entrez_id = ENTREZID)

signatures_data <- bind_rows(
    filter(signatures_data, !is.na(gene_entrez_id)),
    signatures_data %>%
        filter(is.na(gene_entrez_id)) %>%
        rows_update(query, by = "gene_symbol")
)

# Overwrite csv file
signatures_data %>% readr::write_csv("data-raw/hacksig_signatures.csv")

# Save R object
usethis::use_data(signatures_data, internal = TRUE, overwrite = TRUE)
