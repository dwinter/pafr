context("file IO")
# Setup

ali_pafr <- read_paf("test_ali.paf")
ali_tib <- read_paf("test_ali.paf", tibble=TRUE)
ali_tagless <- read_paf("tagless.paf")


tag_names <- c("NM", "ms", "AS", "nn", "tp", "cm", "s1", "s2", "dv", "cg", "zd")
tag_classes <- c("numeric", "numeric", "numeric", "numeric", "character", "numeric", "numeric", "numeric", "numeric", "character", "numeric")

test_that("paf files can be read", {
    #should be a pafr object, with inherits from data.frame
    expect_that(ali_pafr, is_a("pafr"))
    expect_that(ali_pafr, is_a("data.frame"))
    # should be a tibble if we ask for one
    expect_that(ali_tib, is_a("tbl_df"))
    
})

test_that("Alignments are properly represented", {
    
    expect_that( nrow(ali_pafr), equals(5))
    expect_that( sum(ali_pafr[["alen"]]), equals(240777))

    expect_that( nrow(ali_tib), equals(5))
    expect_that( sum(ali_tib[["alen"]]), equals(240777))
})

test_that("alignments include tags", {
    expect_that( names(ali_pafr)[13:23], equals(tag_names))
    expect_that( names(ali_tib)[13:23], equals(tag_names))
})

test_that("tags have right class", {
    expect_that(unname(sapply(ali_pafr, class)[13:23]), equals(tag_classes))
    expect_that(unname(sapply(ali_tib, class)[13:23]), equals(tag_classes))
})

test_that("pafr print fxns work", {
    expect_output(print(ali_pafr), "pafr object with 5 alignments")
    expect_output(print(ali_pafr), "11 tags")
    expect_output(print(ali_pafr), "4 target seqs")
})

test_that("single ali print fxn", {
    expect_output(print(ali_pafr[1,], "Single pafr alignment:\n Q_chr1:660794-665479 v T_chr2:5224220-5228905"))
})

test_that("can handle tagless alignments", {
    expect_equal(ncol(ali_tagless), 12)
    expect_output(print(ali_tagless, "0 tags"))
})

