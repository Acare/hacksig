
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hacksig

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/hacksig)](https://CRAN.R-project.org/package=hacksig)
[![Codecov test
coverage](https://codecov.io/gh/Acare/hacksig/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Acare/hacksig?branch=master)
<!-- badges: end -->

The goal of `hacksig` is to provide a simple and tidy interface to
compute gene expression signature scores in a number of ways. Several
gene signatures and methods are implemented, including:

-   [CINSARC](https://doi.org/10.1038/nm.2174)
-   [ESTIMATE](https://doi.org/10.1038/ncomms3612)
-   [Immunophenoscore](https://doi.org/10.1016/j.celrep.2016.12.019)
-   [Cytolitic activity](https://doi.org/10.1016/j.cell.2014.12.033)
-   and many more (see `?signatures_data` for a complete list of gene
    signatures implemented)

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
library(tidyr)
library(purrr)
library(tibble)
library(future)
library(msigdbr)
```

### Check your signatures

The first thing you should do before computing signature scores is to
check how much of the genes composing a signature are present in your
gene expression matrix. To accomplish this, we can use `check_sig()`:

``` r
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
#> # A tibble: 31 × 4
#>    signature_id     n_genes n_present frac_present
#>    <chr>              <int>     <int>        <dbl>
#>  1 ips_b2m                1         1        1    
#>  2 ips_pd_1               1         1        1    
#>  3 ips_tigit              1         1        1    
#>  4 ips_pd_l2              1         1        1    
#>  5 ips_cd27               1         1        1    
#>  6 ips_icos               1         1        1    
#>  7 ips_act_cd8           26        17        0.654
#>  8 ips_mdsc              20        13        0.65 
#>  9 she2020_irgs          27        17        0.630
#> 10 ayers2017_immexp      18        10        0.556
#> # … with 21 more rows
```

The set of available keywords associated to gene signatures can be
retrieved with:

``` r
get_sig_keywords()
#>  [1] "ayers2017_immexp"   "checkpoint"         "cinsarc"           
#>  [4] "cytolitic activity" "estimate"           "fang2020_irgs"     
#>  [7] "ifng"               "immune"             "immunophenoscore"  
#> [10] "interferon gamma"   "ips_act_cd4"        "ips_act_cd8"       
#> [13] "ips_b2m"            "ips_cd27"           "ips_ctla_4"        
#> [16] "ips_hla_a"          "ips_hla_b"          "ips_hla_dpa1"      
#> [19] "ips_hla_dpb1"       "ips_hla_e"          "ips_hla_f"         
#> [22] "ips_icos"           "ips_lag3"           "ips_mdsc"          
#> [25] "ips_pd_1"           "ips_pd_l1"          "ips_pd_l2"         
#> [28] "ips_tap1"           "ips_tap2"           "ips_tem_cd4"       
#> [31] "ips_tem_cd8"        "ips_tigit"          "ips_tim3"          
#> [34] "ips_treg"           "metabolism"         "metastasis"        
#> [37] "muro2016_ifng "     "rooney2015_cytoact" "sarcoma"           
#> [40] "she2020_irgs"       "stromal"            "sts"               
#> [43] "wu2020_metabolic"
```

We can also check for signatures not implemented in hacksig, that is
custom signatures. For example, we can use the `msigdbr` package to
download the *Hallmark* gene set collection.

``` r
hallmark_list <- msigdbr(species = "Homo sapiens", category = "H") %>%
    distinct(gs_name, gene_symbol) %>%
    nest(genes = c(gene_symbol)) %>%
    mutate(genes = map(genes, compose(as_vector, unname))) %>%
    deframe()
check_sig(test_expr, hallmark_list)
#> # A tibble: 50 × 4
#>    signature_id                        n_genes n_present frac_present
#>    <chr>                                 <int>     <int>        <dbl>
#>  1 HALLMARK_WNT_BETA_CATENIN_SIGNALING      42        27        0.643
#>  2 HALLMARK_APICAL_SURFACE                  44        28        0.636
#>  3 HALLMARK_BILE_ACID_METABOLISM           112        70        0.625
#>  4 HALLMARK_NOTCH_SIGNALING                 32        20        0.625
#>  5 HALLMARK_PI3K_AKT_MTOR_SIGNALING        105        65        0.619
#>  6 HALLMARK_PEROXISOME                     104        64        0.615
#>  7 HALLMARK_COMPLEMENT                     200       123        0.615
#>  8 HALLMARK_TGF_BETA_SIGNALING              54        33        0.611
#>  9 HALLMARK_ESTROGEN_RESPONSE_LATE         200       122        0.61 
#> 10 HALLMARK_INTERFERON_ALPHA_RESPONSE       97        59        0.608
#> # … with 40 more rows
```

### Compute single sample scores

The main function of the package is `hack_sig()` and it permits to
obtain single sample scores from gene signatures. By default,
`hack_sig()` will compute scores for all the signatures implemented in
the package with the original publication method (e.g. weighted sum of
expression values and model coefficients).

``` r
hack_sig(test_expr)
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
#> # A tibble: 20 × 6
#>    sample_id ayers2017_immexp fang2020_irgs muro2016_ifng she2020_irgs
#>    <chr>                <dbl>         <dbl>         <dbl>        <dbl>
#>  1 sample1               5.71         2.96           7.30      1.27   
#>  2 sample2               5.88         2.54           8.36     -1.50   
#>  3 sample3               8.12         6.82           8.87      0.00172
#>  4 sample4               7.83         2.77           8.14     -1.44   
#>  5 sample5               7.89         7.65           8.90      0.165  
#>  6 sample6               9.59         7.83          11.8      -0.573  
#>  7 sample7               8.23         2.24           6.58     -0.00830
#>  8 sample8               7.60         5.46           6.69     -1.40   
#>  9 sample9               7.71         2.02           9.50     -1.81   
#> 10 sample10              6.63         3.41           5.28      0.322  
#> 11 sample11              8.05         4.02           6.22     -0.320  
#> 12 sample12              8.56         8.45           9.02      0.357  
#> 13 sample13              6.35         2.00           7.23      0.453  
#> 14 sample14              8.13         0.988         10.6      -1.65   
#> 15 sample15              7.02         4.09           9.58     -0.990  
#> 16 sample16              7.04         2.31           5.76     -1.08   
#> 17 sample17              7.25         1.49           6.52     -1.04   
#> 18 sample18              7.78         4.21           8.11     -0.453  
#> 19 sample19              9.36         4.20           8.87     -0.320  
#> 20 sample20              8.05        10.5            7.12      3.33   
#> # … with 1 more variable: wu2020_metabolic <dbl>
```

You can also filter for specific signatures (e.g. the immune and stromal
ESTIMATE signatures) and choose a particular single sample method:

``` r
hack_sig(test_expr, signatures = "estimate", method = "zscore")
#> # A tibble: 20 × 3
#>    sample_id estimate_immune estimate_stromal
#>    <chr>               <dbl>            <dbl>
#>  1 sample1           -2.61            -0.240 
#>  2 sample2            1.19            -0.716 
#>  3 sample3           -0.456            0.260 
#>  4 sample4           -0.706            2.18  
#>  5 sample5           -1.06             0.132 
#>  6 sample6           -0.603           -0.0902
#>  7 sample7            0.374           -1.53  
#>  8 sample8           -1.54             0.640 
#>  9 sample9           -1.72            -0.206 
#> 10 sample10           1.32             0.324 
#> 11 sample11           1.49            -0.976 
#> 12 sample12           1.64            -1.20  
#> 13 sample13          -0.530           -0.764 
#> 14 sample14          -0.382            0.386 
#> 15 sample15           1.94            -1.39  
#> 16 sample16          -0.445            0.618 
#> 17 sample17          -0.0235          -1.30  
#> 18 sample18           0.267            1.21  
#> 19 sample19           0.889            0.318 
#> 20 sample20           0.970            2.35
```

Valid choices for single sample `method`s are:

-   `zscore`, for the combined z-score;
-   `ssgsea`, for the single sample GSEA;
-   `singscore`, for the singscore method.

Run `?hack_sig` to see references for these methods.

As in `check_sig()`, the argument `signatures` can also be a list of
gene signatures. For example, we can compute normalized single sample
GSEA scores for the Hallmark gene sets:

``` r
hack_sig(test_expr, hallmark_list, 
         method = "ssgsea", sample_norm = "separate", alpha = 0.5)
#> # A tibble: 20 × 51
#>    sample_id HALLMARK_ADIPOG… HALLMARK_ALLOGR… HALLMARK_ANDROG… HALLMARK_ANGIOG…
#>    <chr>                <dbl>            <dbl>            <dbl>            <dbl>
#>  1 sample1              0.689            0.419            0.946           0.445 
#>  2 sample2              0.500            0.554            0.660           0.616 
#>  3 sample3              0.388            0.367            0.694           0.397 
#>  4 sample4              0.709            0.522            0.797           0.503 
#>  5 sample5              0.340            0.431            0.924           0.0575
#>  6 sample6              0.412            0.531            0.549           0.581 
#>  7 sample7              0.576            0.519            0.891           0.352 
#>  8 sample8              0.488            0.409            0.231           0.250 
#>  9 sample9              0.391            0.300            0.733           0.642 
#> 10 sample10             0.381            0.805            0.308           0.292 
#> 11 sample11             0.248            0.757            0.650           1     
#> 12 sample12             0.999            1                0.971           0.123 
#> 13 sample13             0.783            0                0.371           0.610 
#> 14 sample14             0.445            0.510            0.754           0.236 
#> 15 sample15             0.519            0.591            0               0.0382
#> 16 sample16             0.450            0.632            0.536           0.0674
#> 17 sample17             1                0.969            1               0     
#> 18 sample18             0.313            0.637            0.542           0.204 
#> 19 sample19             0.587            0.126            0.559           0.289 
#> 20 sample20             0                0.654            0.724           0.382 
#> # … with 46 more variables: HALLMARK_APICAL_JUNCTION <dbl>,
#> #   HALLMARK_APICAL_SURFACE <dbl>, HALLMARK_APOPTOSIS <dbl>,
#> #   HALLMARK_BILE_ACID_METABOLISM <dbl>,
#> #   HALLMARK_CHOLESTEROL_HOMEOSTASIS <dbl>, HALLMARK_COAGULATION <dbl>,
#> #   HALLMARK_COMPLEMENT <dbl>, HALLMARK_DNA_REPAIR <dbl>,
#> #   HALLMARK_E2F_TARGETS <dbl>,
#> #   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <dbl>, …
```

There are three methods for which `hack_sig()` cannot be used to compute
gene signature scores with the original method. These are: CINSARC,
ESTIMATE and the Immunophenoscore.

For the CINSARC classification, you must provide a vector with distant
metastasis status:

``` r
rand_dm <- sample(c(0, 1), size = ncol(test_expr), replace = TRUE)
hack_cinsarc(test_expr, rand_dm)
#> # A tibble: 20 × 2
#>    sample_id cinsarc_class
#>    <chr>     <chr>        
#>  1 sample1   C2           
#>  2 sample2   C2           
#>  3 sample3   C2           
#>  4 sample4   C2           
#>  5 sample5   C2           
#>  6 sample6   C1           
#>  7 sample7   C1           
#>  8 sample8   C1           
#>  9 sample9   C1           
#> 10 sample10  C1           
#> 11 sample11  C1           
#> 12 sample12  C2           
#> 13 sample13  C2           
#> 14 sample14  C1           
#> 15 sample15  C1           
#> 16 sample16  C1           
#> 17 sample17  C2           
#> 18 sample18  C1           
#> 19 sample19  C2           
#> 20 sample20  C1
```

Immune, stromal, ESTIMATE and tumor purity scores from the ESTIMATE
method can be obtained with:

``` r
hack_estimate(test_expr)
#> # A tibble: 20 × 5
#>    sample_id immune_score stroma_score estimate_score purity_score
#>    <chr>            <dbl>        <dbl>          <dbl>        <dbl>
#>  1 sample1          -634.         775.           141.        0.811
#>  2 sample2          2121.         703.          2824.        0.524
#>  3 sample3           727.         800.          1526.        0.676
#>  4 sample4           714.        2029.          2744.        0.534
#>  5 sample5           183.        1121.          1304.        0.699
#>  6 sample6          1192.        1166.          2358.        0.581
#>  7 sample7          1323.         377.          1699.        0.657
#>  8 sample8           546.        1155.          1701.        0.656
#>  9 sample9           323.        1156.          1479.        0.681
#> 10 sample10         1620.        1300.          2920.        0.512
#> 11 sample11         2044.         530.          2574.        0.555
#> 12 sample12         1836.         777.          2613.        0.550
#> 13 sample13          622.         783.          1405.        0.689
#> 14 sample14         1171.        1003.          2174.        0.603
#> 15 sample15         2396.         434.          2829.        0.523
#> 16 sample16         1303.        1274.          2578.        0.554
#> 17 sample17         1189.         691.          1880.        0.636
#> 18 sample18          857.        1515.          2372.        0.579
#> 19 sample19         1658.         982.          2640.        0.547
#> 20 sample20         1642.        2487.          4129.        0.352
```

Finally, the raw immunophenoscore and its discrete counterpart can be
obtained with:

``` r
hack_immunophenoscore(test_expr)
#> # A tibble: 20 × 3
#>    sample_id raw_score ips_score
#>    <chr>         <dbl>     <dbl>
#>  1 sample1      0.930          3
#>  2 sample2     -0.224          0
#>  3 sample3      0.0827         0
#>  4 sample4     -0.334          0
#>  5 sample5      1.63           5
#>  6 sample6      1.71           6
#>  7 sample7     -0.148          0
#>  8 sample8      0.506          2
#>  9 sample9     -0.233          0
#> 10 sample10    -2.32           0
#> 11 sample11     0.112          0
#> 12 sample12     1.44           5
#> 13 sample13    -0.863          0
#> 14 sample14     0.907          3
#> 15 sample15     0.637          2
#> 16 sample16     0.746          2
#> 17 sample17     0.334          1
#> 18 sample18    -0.322          0
#> 19 sample19    -1.09           0
#> 20 sample20     0.806          3
```

### Classify your samples

If you want to categorize your samples into two or more signature
classes based on a score cutoff, you can use `hack_class()` after
`hack_sig()`:

``` r
test_expr %>% 
    hack_sig("estimate", method = "singscore", direction = "up") %>% 
    hack_class()
#> # A tibble: 20 × 3
#>    sample_id estimate_immune estimate_stromal
#>    <chr>     <chr>           <chr>           
#>  1 sample1   low             low             
#>  2 sample2   high            low             
#>  3 sample3   low             low             
#>  4 sample4   low             high            
#>  5 sample5   low             high            
#>  6 sample6   low             high            
#>  7 sample7   high            low             
#>  8 sample8   low             high            
#>  9 sample9   low             low             
#> 10 sample10  high            high            
#> 11 sample11  high            low             
#> 12 sample12  high            low             
#> 13 sample13  low             low             
#> 14 sample14  high            high            
#> 15 sample15  high            low             
#> 16 sample16  low             high            
#> 17 sample17  high            low             
#> 18 sample18  low             high            
#> 19 sample19  high            high            
#> 20 sample20  high            high
```

### Speed-up computation time

`hacksig` supports the `future` framework for parallel computation. Our
single sample method implementations based on ranks (i.e. single sample
GSEA and singscore) are slower than their counterparts implemented in
`GSVA` and `singscore`. Hence, to speed-up computation time you can use
the `future` package:

``` r
plan(multisession)
hack_sig(test_expr, hallmark_list, method = "ssgsea")
#> # A tibble: 20 × 51
#>    sample_id HALLMARK_ADIPOG… HALLMARK_ALLOGR… HALLMARK_ANDROG… HALLMARK_ANGIOG…
#>    <chr>                <dbl>            <dbl>            <dbl>            <dbl>
#>  1 sample1              1624.             711.           2005.            1376. 
#>  2 sample2              1280.            1052.           1125.            1673. 
#>  3 sample3              1066.             666.           1129.            1080. 
#>  4 sample4              1614.             954.           1679.            1575. 
#>  5 sample5               826.             835.           1877.             610. 
#>  6 sample6              1025.            1045.            619.            1514. 
#>  7 sample7              1436.            1020.           1787.            1273. 
#>  8 sample8              1157.             938.           -103.             760. 
#>  9 sample9              1006.             499.           1438.            2046. 
#> 10 sample10              735.            1514.             66.9            873. 
#> 11 sample11              574.            1500.           1091.            2737. 
#> 12 sample12             2432.            1959.           2207.              20.1
#> 13 sample13             1824.             113.            163.            1659. 
#> 14 sample14             1179.            1059.           1424.             755. 
#> 15 sample15             1060.            1043.           -930.             130. 
#> 16 sample16             1098.            1208.            601.             527. 
#> 17 sample17             2396.            1900.           2058.             397. 
#> 18 sample18              808.            1225.            837.             924. 
#> 19 sample19             1167.             300.            814.            1275. 
#> 20 sample20              145.            1262.           1318.             859. 
#> # … with 46 more variables: HALLMARK_APICAL_JUNCTION <dbl>,
#> #   HALLMARK_APICAL_SURFACE <dbl>, HALLMARK_APOPTOSIS <dbl>,
#> #   HALLMARK_BILE_ACID_METABOLISM <dbl>,
#> #   HALLMARK_CHOLESTEROL_HOMEOSTASIS <dbl>, HALLMARK_COAGULATION <dbl>,
#> #   HALLMARK_COMPLEMENT <dbl>, HALLMARK_DNA_REPAIR <dbl>,
#> #   HALLMARK_E2F_TARGETS <dbl>,
#> #   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <dbl>, …
```

## Contributing

If you have any suggestions about adding new features to `hacksig`,
please open an issue request on
[GitHub](https://github.com/Acare/hacksig/issues). Gene-level
information about gene signatures are stored in the R object
`signatures_data` and can be used as a template for requests.
