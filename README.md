
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hacksig

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/hacksig)](https://CRAN.R-project.org/package=hacksig)
[![Codecov test
coverage](https://codecov.io/gh/Acare/hacksig/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Acare/hacksig?branch=master)
<!-- badges: end -->

The goal of `hacksig` is to provide a simple and tidy interface to
compute single sample scores for a manually curated collection of gene
signatures and methods applied in cancer transcriptomics.

Several gene signatures and methods are implemented, including:

-   [CINSARC](https://doi.org/10.1038/nm.2174);
-   [ESTIMATE](https://doi.org/10.1038/ncomms3612);
-   [Immunophenoscore](https://doi.org/10.1016/j.celrep.2016.12.019);
-   [Cytolitic activity](https://doi.org/10.1016/j.cell.2014.12.033);
-   and many more (use `get_sig_info()` for a complete list of gene
    signatures implemented)

One can choose to apply either the original publication method or one of
three single sample scoring alternatives, namely: combined z-score,
single sample GSEA and singscore.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Acare/hacksig")
```

## Usage

To show example usage of the package, we will use a simulated gene
expression matrix `test_expr`, which is stored as an R object in
`hacksig`.

``` r
library(hacksig)
library(dplyr)
library(future)
```

### Get available signatures

``` r
get_sig_info()
#> # A tibble: 35 × 4
#>   signature_id     signature_keywords             publication_doi    description
#>   <chr>            <chr>                          <chr>              <chr>      
#> 1 ayers2017_immexp ayers2017_immexp|immune        10.1172/JCI91190   <NA>       
#> 2 cinsarc          cinsarc|sarcoma|sts|metastasis 10.1038/nm.2174    <NA>       
#> 3 estimate_immune  estimate|immune                10.1038/ncomms3612 <NA>       
#> 4 estimate_stromal estimate|stromal               10.1038/ncomms3612 <NA>       
#> 5 fang2020_irgs    fang2020_irgs|immune           10.18632/aging.20… <NA>       
#> # … with 30 more rows
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

### Classify your samples

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
#> x rooney2015_cytoact
#> x ips_tap1
#> x ips_tap2
#> x ips_hla_a
#> x ips_hla_b
#> x ips_hla_c
#> x ips_hla_dpa1
#> x ips_hla_dpb1
#> x ips_hla_e
#> x ips_hla_f
#> x ips_ctla_4
#> x ips_lag3
#> x ips_tim3
#> x ips_pd_l1
#> x ips_ido1
#> # A tibble: 20 × 15
#>   sample_id ayers2017_immexp cinsarc estimate_immune estimate_stromal
#>   <chr>                <dbl>   <dbl>           <dbl>            <dbl>
#> 1 sample1             -3914.   -13.5           -636.             778.
#> 2 sample2             -3348. -1070.            2118.             703.
#> 3 sample3              1697.  1805.             725.             805.
#> 4 sample4               366.   326.             737.            2031.
#> 5 sample5               969.   290.             181.            1129.
#> # … with 15 more rows, and 10 more variables: fang2020_irgs <dbl>,
#> #   ips_act_cd4 <dbl>, ips_act_cd8 <dbl>, ips_mdsc <dbl>, ips_tem_cd4 <dbl>,
#> #   ips_tem_cd8 <dbl>, ips_treg <dbl>, muro2016_ifng <dbl>, she2020_irgs <dbl>,
#> #   wu2020_metabolic <dbl>
```

## Contributing

If you have any suggestions about adding new features to `hacksig`,
please open an issue request on
[GitHub](https://github.com/Acare/hacksig/issues). Gene-level
information about gene signatures are stored in
`data-raw/hacksig_signatures.csv` and can be used as a template for
requests.
