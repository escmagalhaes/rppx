
# rppx

<!-- badges: start -->
<!-- badges: end -->

Designed to work with pre-normalized data derived from Reverse Phase Protein 
Arrays (RPPA), downstream of RPPASPACE. Provides a compendium of functions 
related to analytical pipelines with graphical visualization and figure generation.

## Installation

You can install the beta version of rppx from [GitHub](https://github.com/) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("escmagalhaes/rppx@v0.1.0-beta")
```
A few dependencies belong to Bioconductor and will be automatically installed  

## Example

This is a basic example which shows you how to use this package:

``` r
library(rppx)
dataframe<-data.frame(
status = c("Alive","alive","Dead ","dead",NA),
grade  = c("low","Low","HIGH","high","medium")
)
get_rawtable_values(data=dataframe,features=c("status","grade"))
```

## Status

This package is under development (beta version)

## License

MIT License - see related LICENSE file
