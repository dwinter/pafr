context("manipulating pafr objects")

ali_notag <- read_paf("tagless.paf", tibble=TRUE)
ali_all_inv <- read_paf("all_inverted.paf", tibble=TRUE)

test_that("chrom-sizes fxn returns data.frame from pafr object", {
    cs <- chrom_sizes(ali_pafr)
    expect_that(cs, is_a("list"))
    expect_length(cs, 2)
})

test_that("chrom-sizes fxn returns data.frame from tibbl object", {
    cs <- chrom_sizes(ali_tib)
    expect_that(cs, is_a("list"))
    expect_length(cs, 2)
})

test_that("filtering by alignment-type works", {
    expect_error(filter_secondary_alignments(ali_notag)) # no tags, no way to filter
    expect_equal(nrow(filter_secondary_alignments(ali_all_inv)), 4) # no tags, no way to filter
    expect_equal(nrow(filter_secondary_alignments(ali_all_inv, remove_inversions=TRUE)), 0) # no tags, no way to filter
})
