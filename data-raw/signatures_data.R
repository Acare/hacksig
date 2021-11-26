library(dplyr)
library(org.Hs.eg.db) # version 3.12.0
library(hgu133plus2.db) # version 3.2.3

# Add some signatures to hacksig_signatures.csv. Then:
signatures_data <- readr::read_csv("data-raw/hacksig_signatures.csv",
                                   col_types = "ccccnccc") %>%
    mutate(
        gene_symbol = case_when(
            gene_symbol == "WISP1" ~ "CCN4",    # estimate
            gene_symbol == "GPR124" ~ "ADGRA2",
            gene_symbol == "LPPR4" ~ "PLPPR4",
            gene_symbol == "TXNDC3" ~ "NME8",
            gene_symbol == "ODZ4" ~ "TENM4",
            gene_symbol == "FYB" ~ "FYB1",
            TRUE ~ gene_symbol
        ),
        gene_symbol = case_when(
            gene_symbol == "IARS" ~ "IARS1",    # immunophenoscore
            gene_symbol == "DARS" ~ "DARS1",
            gene_symbol == "SDPR" ~ "CAVIN2",
            TRUE ~ gene_symbol
        ))

# Add entrez ids from gene symbols
query <- select(org.Hs.eg.db,
                keys = unique(signatures_data$gene_symbol),
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL")

# add gene symbols and entrez ids from probe ids
query <- select(hgu133plus2.db,
                keys = "210220_at",
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "PROBEID")

# signatures_data <- signatures_data %>%
#     left_join(query, by = c("gene_symbol" = "SYMBOL")) %>%
#     mutate(entrez_gene_id = ENTREZID, .keep = "unused")

signatures_data %>% readr::write_csv("data-raw/hacksig_signatures.csv")

usethis::use_data(signatures_data, internal = TRUE, overwrite = TRUE)
