context("Test the coverage plots code")

ali_pafr <- read_paf("test_ali.paf")


test_that("plot_coverage produces a plot", {
    p <- plot_coverage(ali)
    expect_is(p, "ggplot")

})

test_that("can shade alignments basde on ali attributes", {
    p_monochrome <- plot_coverage(ali)
    p_technocolour <- plot_coverage(ali, fill="qname") 
    expect_true('fill' %in% names(p_technocolour$layers[[2]]$mapping))
    expect_false('fill' %in% names(p_monochrome$layers[[2]]$mapping))
})

test_that("We can plot query sequences", {
    p_q <- plot_coverage(ali, target=FALSE)
    expect_is(p_q, "ggplot")
    expect_equal(layer_data(p_q)[1,2], 7744434)
})
