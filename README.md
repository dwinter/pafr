---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# pafr

Read, manipulate and visualize 'Pairwise mApping Format' data in R

## Under contstruction

This package is in the process of being turned from reseach code into a nice
usable package. Much of the code is already here, but documentation tests and
added-flexibility are still being added. 

Nevertheless, you can check out the code and run it for youself. Here's how to
make a basic whole-genome dotplot from a PAF file:

```r
library(pafr)
ali <- read_paf("my_paf_file.paf")
dotplot(ali)
```
