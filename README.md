<!-- README.md is generated from README.Rmd. Please edit that file -->
![Devlopment stage](https://img.shields.io/badge/development--stage-alpha-lightgrey) [![Travis-CI Build Status](https://travis-ci.org/dwinter/pafr.svg?branch=master)](https://travis-ci.org/dwinter/pafr) [![Coverage Status](https://img.shields.io/codecov/c/github/dwinter/pafr/master.svg)](https://codecov.io/github/dwinter/pafr?branch=master)

pafr
====

Read, manipulate and visualize 'Pairwise mApping Format' data in R

Under contstruction
-------------------

This package is in the process of being turned from reseach code into a nice usable package. Much of the code is already here, but documentation tests and added-flexibility are still being added.

Nevertheless, you can check out the code and run it for youself. Here's how to make a basic whole-genome dotplot from a PAF file:

``` r
library(pafr)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following object is masked from 'package:MoreUtils':
#> 
#>     last
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> Loading required package: tibble
test_alignment <- system.file("extdata", "fungi.paf", package="pafr")
ali <- read_paf(test_alignment)
dotplot(ali)
```

![](man/figures/README-dotplot-1.png)

Installation
------------

If you want to play with the code already here, you can install this package using devtools

``` r
#install.packages(pafr)
devtools::install_github("dwinter/pafr")
```

Bugs/Issues
-----------

Please use the issue tracker on this repo to let us know about any bugs or issues that arise as you use this package
