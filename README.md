
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteins

<!-- badges: start -->

<!-- badges: end -->

This package contains the `proteins` dataset and some utility functions
related to protein analyses. It accompanies a workshop-style class that
provides an introduction to the emerging field of Data Science in R,
including data analysis and visualization, with a particular focus on
its utility for biological insight.

## Installation

You **cannot** yet install the released version of tidybiology from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("proteins")
```

So in the meantime, use the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hirscheylab/proteins")
```

## Example

This is how to take a `glimpse` into the proteins dataset:

``` r
library(proteins)
glimpse(proteins)
#> Observations: 20,430
#> Variables: 8
#> $ uniprot_id       <chr> "P04217", "Q9NQ94", "P01023", "A8K2U0", "U3KPV4…
#> $ gene_name        <chr> "A1BG", "A1CF", "A2M", "A2ML1", "A3GALT2", "A4G…
#> $ gene_name_alt    <chr> NA, "ACF ASP", "CPAMD5 FWP007", "CPAMD9", "A3GA…
#> $ protein_name     <chr> "Alpha-1B-glycoprotein ", "APOBEC1 complementat…
#> $ protein_name_alt <chr> "Alpha-1-B glycoprotein)", "APOBEC1-stimulating…
#> $ sequence         <chr> "MSMLVVFLLLWGVTWGPVTEAAIFYETQPSLWAESESLLKPLANVT…
#> $ length           <dbl> 495, 594, 1474, 1454, 340, 353, 340, 546, 672, …
#> $ mass             <dbl> 54254, 65202, 163291, 161107, 38754, 40499, 394…
```
