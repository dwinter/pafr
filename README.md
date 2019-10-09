<!-- README.md is generated from README.Rmd. Please edit that file -->
![Devlopment stage](https://img.shields.io/badge/development--stage-beta-orange) [![Travis-CI Build Status](https://travis-ci.org/dwinter/pafr.svg?branch=master)](https://travis-ci.org/dwinter/pafr) [![Coverage Status](https://img.shields.io/codecov/c/github/dwinter/pafr/master.svg)](https://codecov.io/github/dwinter/pafr?branch=master)

pafr
====

Read, manipulate and visualize 'Pairwise mApping Format' data in R

Under contstruction
-------------------

This package is in the process of being turned from reseach code into a nice usable package. The package is now feature-complete, but is still no entirely documented and likely has bugs that need to be fixed.

Nevertheless, you can check out the code and run it for youself. Here's how to make a basic whole-genome dotplot from a PAF file:

``` r
library(pafr, quietly=TRUE)
test_alignment <- system.file("extdata", "fungi.paf", package="pafr")
ali <- read_paf(test_alignment)
dotplot(ali, label_seqs=TRUE)
```

![](man/figures/README-dotplot-1.png)

Installation
------------

If you want to play with the code already here, you can install this package using devtools

``` r
#install.packages(devtools)
devtools::install_github("dwinter/pafr")
```

Bugs/Issues
-----------

Please use the issue tracker on this repo to let us know about any bugs or issues that arise as you use this package.
