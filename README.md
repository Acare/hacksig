
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hacksig

<!-- badges: start -->
<!-- badges: end -->

The goal of hacksig is to …

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Acare/hacksig")
```

## Getting started

The first thing one should do before computing signature scores is to
check how much of the genes composing a signature are present in your
gene expression matrix. To accomplish this, we can use `check_sig()`:

``` r
library(hacksig)
check_sig(test_expr)
#> # A tibble: 35 × 4
#>    signature_id     n_genes n_present frac_present
#>    <chr>              <int>     <int>        <dbl>
#>  1 ips_b2m                1         1        1    
#>  2 ips_pd_1               1         1        1    
#>  3 ips_tigit              1         1        1    
#>  4 ips_pd_l2              1         1        1    
#>  5 ips_cd27               1         1        1    
#>  6 ips_icos               1         1        1    
#>  7 muro2016_ifng          6         4        0.667
#>  8 wu2020_metabolic      30        20        0.667
#>  9 ips_act_cd8           26        17        0.654
#> 10 ips_mdsc              20        13        0.65 
#> # … with 25 more rows
```

By default, `check_sig()` will compute statistics for every signature
implemented in `hacksig`. You can filter for specific signatures by
entering a keyword in the `signatures` argument:

``` r
check_sig(test_expr, signatures = "immune")
#> # A tibble: 1 × 4
#>   signature_id    n_genes n_present frac_present
#>   <chr>             <int>     <int>        <dbl>
#> 1 estimate_immune     141        74        0.525
```

The set of available keywords associated to gene signatures can be
retrieved with:

``` r
get_sig_keywords()
#>  [1] "act cd4"            "act cd8"            "b2m"               
#>  [4] "cd27"               "checkpoint"         "ctla-4"            
#>  [7] "ctla4"              "cytolitic activity" "estimate"          
#> [10] "hla-a"              "hla-b"              "hla-c"             
#> [13] "hla-dpa1"           "hla-dpb1"           "hla-e"             
#> [16] "hla-f"              "ici"                "icos"              
#> [19] "ido1"               "ifng"               "immune"            
#> [22] "immunophenoscore"   "interferon gamma"   "ips"               
#> [25] "lag3"               "mdsc"               "metabolism"        
#> [28] "metastasis"         "pd-1"               "pd-l1"             
#> [31] "pd-l2"              "pd1"                "pdl1"              
#> [34] "pdl2"               "sarcoma"            "stromal"           
#> [37] "sts"                "tap1"               "tap2"              
#> [40] "tem cd4"            "tem cd8"            "tigit"             
#> [43] "tim3"               "treg"
```
