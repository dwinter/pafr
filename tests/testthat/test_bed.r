context("Test bed file reading")

test_that("Bed files can be read", {
    B <- read_bed("test_I.bed")
    expect_equal(ncol(B), 3)
    expect_equal(names(B), c("chrom", "start", "end"))
    expect_equal(B[["chrom"]], c("Q_chr1", "Q_chr2"))
})

test_that("Error on bad bed", {
    expect_error(B <- read_bed("bad_I.bed"))
})
