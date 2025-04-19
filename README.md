
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hacksig <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/hacksig)](https://CRAN.R-project.org/package=hacksig)
[![Codecov test
coverage](https://codecov.io/gh/Acare/hacksig/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Acare/hacksig?branch=master)
[![R-CMD-check](https://github.com/Acare/hacksig/workflows/R-CMD-check/badge.svg)](https://github.com/Acare/hacksig/actions)
[![R-CMD-check](https://github.com/Acare/hacksig/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Acare/hacksig/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of `hacksig` is to provide a simple and tidy interface to
compute single sample scores for gene signatures and methods applied in
cancer transcriptomics.

Scores can be obtained either for custom lists of genes or for a
manually curated collection of gene signatures, including:

- [CINSARC](https://doi.org/10.1038/nm.2174);
- [ESTIMATE](https://doi.org/10.1038/ncomms3612);
- [Immunophenoscore](https://doi.org/10.1016/j.celrep.2016.12.019);
- [Cytolitic activity](https://doi.org/10.1016/j.cell.2014.12.033);
- and more (use `get_sig_info()` to get a complete list of the
  implemented signatures)

At present, signature scores can be obtained either with the original
publication method or using one of three single sample scoring
alternatives, namely: *combined z-score*, *single sample GSEA* and
*singscore*.

## Installation

You can install the last stable version of hacksig from
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

### Available signatures

``` r
get_sig_info()
#> # A tibble: 40 × 4
#>   signature_id       signature_keywords              publication_doi description
#>   <chr>              <chr>                           <chr>           <chr>      
#> 1 ayers2017_immexp   ayers2017_immexp|immune expand… 10.1172/JCI911… Immune exp…
#> 2 bai2019_immune     bai2019_immune|head and neck s… 10.1155/2019/3… Immune/inf…
#> 3 cinsarc            cinsarc|metastasis|sarcoma|sts  10.1038/nm.2174 Biomarker …
#> 4 dececco2014_int172 dececco2014_int172|head and ne… 10.1093/annonc… Signature …
#> 5 eschrich2009_rsi   eschrich2009_rsi|radioresistan… 10.1016/j.ijro… Genes aime…
#> # ℹ 35 more rows
```

### Check your signatures

``` r
check_sig(test_expr, signatures = "estimate")
#> # A tibble: 2 × 5
#>   signature_id     n_genes n_present frac_present missing_genes
#>   <chr>              <int>     <int>        <dbl> <list>       
#> 1 estimate_stromal     141        91        0.645 <chr [50]>   
#> 2 estimate_immune      141        74        0.525 <chr [67]>
```

### Compute single sample scores

``` r
hack_sig(test_expr, signatures = c("ifng", "cinsarc"), method = "zscore")
#> # A tibble: 20 × 3
#>   sample_id cinsarc muro2016_ifng
#>   <chr>       <dbl>         <dbl>
#> 1 sample1   -0.482         -0.511
#> 2 sample10  -0.0926        -1.60 
#> 3 sample11   0.730         -1.03 
#> 4 sample12  -0.625          0.851
#> 5 sample13   0.930         -0.369
#> # ℹ 15 more rows
```

### Stratify your samples

``` r
test_expr %>% 
    hack_sig("estimate", method = "singscore", direction = "up") %>% 
    stratify_sig(cutoff = "median")
#> # A tibble: 20 × 3
#>   sample_id estimate_immune estimate_stromal
#>   <chr>     <chr>           <chr>           
#> 1 sample1   low             low             
#> 2 sample10  high            high            
#> 3 sample11  high            low             
#> 4 sample12  high            low             
#> 5 sample13  low             low             
#> # ℹ 15 more rows
```

### Speed-up computation time

``` r
plan(multisession)
hack_sig(test_expr, method = "ssgsea")
#> Warning: ℹ No genes are present in 'expr_data' for the following signatures:
#> ✖ zhu2021_ferroptosis
#> ✖ rooney2015_cyt
#> # A tibble: 20 × 39
#>   sample_id ayers2017_immexp bai2019_immune cinsarc dececco2014_int172
#>   <chr>                <dbl>          <dbl>   <dbl>              <dbl>
#> 1 sample1             -3914.          2316.   -13.5              1288.
#> 2 sample10             1077.           575.   801.                811.
#> 3 sample11              501.          -490.  1340.               1244.
#> 4 sample12             2315.          1034.  -151.                981.
#> 5 sample13            -2179.           327.  1737.               1288.
#> # ℹ 15 more rows
#> # ℹ 34 more variables: eschrich2009_rsi <dbl>, estimate_immune <dbl>,
#> #   estimate_stromal <dbl>, eustace2013_hypoxia <dbl>,
#> #   fan2021_ferroptosis <dbl>, fang2021_irgs <dbl>, han2021_ferroptosis <dbl>,
#> #   he2021_ferroptosis_a <dbl>, he2021_ferroptosis_b <dbl>, hu2021_derbp <dbl>,
#> #   huang2022_ferroptosis <dbl>, ips_cp <dbl>, ips_ec <dbl>, ips_mhc <dbl>,
#> #   ips_sc <dbl>, li2021_ferroptosis_a <dbl>, li2021_ferroptosis_b <dbl>, …
```

## Contributing

If you have any suggestions about adding new features or signatures to
`hacksig`, please create an issue on
[GitHub](https://github.com/Acare/hacksig/issues). Gene-level
information about gene signatures are stored in
`data-raw/hacksig_signatures.csv` and can be used as a template for
requests.
