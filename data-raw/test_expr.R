library(dplyr)
library(org.Hs.eg.db) # version 3.12.0

filt_gene_symbols <- keys(org.Hs.eg.db, "SYMBOL") %>%
    as_tibble() %>%
    filter(!stringr::str_detect(value, "^LINC"),
           !stringr::str_detect(value, "^LOC"),
           !stringr::str_detect(value, "-")) %>%
    pull(value)

set.seed(123)

test_expr <- matrix(rnorm(20000 * 20, 8, 3), nrow = 20000, ncol = 20)
test_expr[test_expr < 0] <- 1
test_expr <- apply(test_expr, MARGIN = 2, FUN = signif, 3)
colnames(test_expr) <- paste0("sample", seq_len(ncol(test_expr)))
rownames(test_expr) <- sample(filt_gene_symbols, size = nrow(test_expr))

usethis::use_data(test_expr, overwrite = TRUE)
