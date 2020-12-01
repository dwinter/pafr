context("file IO")
# Setup

ali_pafr <- read_paf("test_ali.paf")
ali_skip_tags <- read_paf("test_ali.paf", include_tags=FALSE)
ali_tib <- read_paf("test_ali.paf", tibble=TRUE)
ali_tagless <- read_paf("tagless.paf")
ali_df <- read.table("tagless.paf")


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
test_that("We can ignore tags read_paf", {
    expect_that(ali_skip_tags, is_a("pafr"))
    expect_that( ncol(ali_skip_tags), equals(12))
})


test_that("can convert data.frame to pafr", {
    expect_that(ali_df, is_a("data.frame"))
    ali_converted <- as_paf(ali_df)
    expect_that(ali_converted, is_a("pafr"))
    
    expect_that( nrow(ali_converted), equals(5))
    expect_that( ncol(ali_converted), equals(12))
    expect_that( sum(ali_converted[["alen"]]), equals(240777))
})

test_that("throw errors on converting data.frames", {
    ali_df$extra_col <- "shouldn't be here"
    expect_error( as_paf(ali_df) )
    ali_df$extra_col <- NULL
    ali_df[[12]] <- "Should be a number"
    expect_error( as_paf(ali_df) )
})
