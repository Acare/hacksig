
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hacksig

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/hacksig)](https://CRAN.R-project.org/package=hacksig)
[![Codecov test
coverage](https://codecov.io/gh/Acare/hacksig/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Acare/hacksig?branch=master)
[![R-CMD-check](https://github.com/Acare/hacksig/workflows/R-CMD-check/badge.svg)](https://github.com/Acare/hacksig/actions)
<!-- badges: end -->

The goal of `hacksig` is to provide a simple and tidy interface to
compute single sample scores for gene signatures and methods applied in
cancer transcriptomics.

Scores can be obtained either for custom lists of genes or for a
manually curated collection of gene signatures, including:

-   [CINSARC](https://doi.org/10.1038/nm.2174);
-   [ESTIMATE](https://doi.org/10.1038/ncomms3612);
-   [Immunophenoscore](https://doi.org/10.1016/j.celrep.2016.12.019);
-   [Cytolitic activity](https://doi.org/10.1016/j.cell.2014.12.033);
-   and more (use `get_sig_info()` for a complete list of gene
    signatures implemented)

One can choose to apply either the original publication method or one of
three single sample scoring alternatives, namely: combined z-score,
single sample GSEA and singscore.

## Installation

You can install the released version of hacksig from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hacksig")
```

Or the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Acare/hacksig")
```

## Usage

You can learn more about usage of the package in `vignette("hacksig")`.

``` r
library(hacksig)
library(dplyr)
library(future)
```

### Get available signatures

``` r
get_sig_info()
#> # A tibble: 23 × 4
#>   signature_id       signature_keywords              publication_doi description
#>   <chr>              <chr>                           <chr>           <chr>      
#> 1 ayers2017_immexp   ayers2017_immexp|immune expand… 10.1172/JCI911… Immune exp…
#> 2 bai2019_immune     bai2019_immune|head and neck|h… 10.1155/2019/3… Immune/inf…
#> 3 cinsarc            cinsarc|metastasis|sarcoma|sts  10.1038/nm.2174 Biomarker …
#> 4 dececco2014_int172 dececco2014_int172|head and ne… 10.1093/annonc… Signature …
#> 5 eschrich2009_rsi   eschrich2009_rsi|radioresistan… 10.1016/j.ijro… Genes aime…
#> # … with 18 more rows
```

### Check your signatures

``` r
check_sig(test_expr, signatures = "estimate")
#> # A tibble: 2 × 5
#>   signature_id     n_genes n_present frac_present missing_genes
#>   <chr>              <int>     <int>        <dbl> <named list> 
#> 1 estimate_stromal     141        91        0.645 <chr [50]>   
#> 2 estimate_immune      141        74        0.525 <chr [67]>
```

### Compute single sample scores

``` r
hack_sig(test_expr, signatures = c("ifng", "cinsarc"), method = "zscore")
#> # A tibble: 20 × 3
#>   sample_id cinsarc muro2016_ifng
#>   <chr>       <dbl>         <dbl>
#> 1 sample1    -0.482       -0.511 
#> 2 sample2    -2.61         0.400 
#> 3 sample3     1.44         0.347 
#> 4 sample4    -0.538        0.0849
#> 5 sample5    -0.537        0.390 
#> # … with 15 more rows
```

### Stratify your samples

``` r
test_expr %>% 
    hack_sig("estimate", method = "singscore", direction = "up") %>% 
    hack_class(cutoff = "median")
#> # A tibble: 20 × 3
#>   sample_id estimate_immune estimate_stromal
#>   <chr>     <chr>           <chr>           
#> 1 sample1   low             low             
#> 2 sample2   high            low             
#> 3 sample3   low             low             
#> 4 sample4   low             high            
#> 5 sample5   low             high            
#> # … with 15 more rows
```

### Speed-up computation time

``` r
plan(multisession)
hack_sig(test_expr, method = "ssgsea")
#> Warning: ℹ No genes are present in 'expr_data' for the following signatures:
#> x rooney2015_cyt
#> # A tibble: 20 × 23
#>   sample_id ayers2017_immexp bai2019_immune cinsarc dececco2014_int172
#>   <chr>                <dbl>          <dbl>   <dbl>              <dbl>
#> 1 sample1             -3914.          2316.   -13.5              1288.
#> 2 sample2             -3348.          1350. -1070.               1322.
#> 3 sample3              1697.          1829.  1805.                685.
#> 4 sample4               366.          5611.   326.               1684.
#> 5 sample5               969.          1224.   290.                718.
#> # … with 15 more rows, and 18 more variables: eschrich2009_rsi <dbl>,
#> #   estimate_immune <dbl>, estimate_stromal <dbl>, eustace2013_hypoxia <dbl>,
#> #   fang2021_irgs <dbl>, hu2021_derbp <dbl>, ips_cp <dbl>, ips_ec <dbl>,
#> #   ips_mhc <dbl>, ips_sc <dbl>, li2021_irgs <dbl>, liu2020_immune <dbl>,
#> #   liu2021_mgs <dbl>, lohavanichbutr2013_hpvneg <dbl>, muro2016_ifng <dbl>,
#> #   qiang2021_irgs <dbl>, she2020_irgs <dbl>, wu2020_metabolic <dbl>
```

## Contributing

If you have any suggestions about adding new features to `hacksig`,
please open an issue request on
[GitHub](https://github.com/Acare/hacksig/issues). Gene-level
information about gene signatures are stored in
`data-raw/hacksig_signatures.csv` and can be used as a template for
requests.
