context("Test the coverage plots code")


test_that("plot_coverage produces a plot", {
    p <- plot_coverage(ali_pafr)
    expect_is(p, "ggplot")

})

test_that("can shade alignments basde on ali attributes", {
    p_monochrome <- plot_coverage(ali_pafr)
    p_technocolour <- plot_coverage(ali_pafr, fill="qname") 
    expect_true('fill' %in% names(p_technocolour$layers[[2]]$mapping))
    expect_false('fill' %in% names(p_monochrome$layers[[2]]$mapping))
})

test_that("We can plot query sequences", {
    p_q <- plot_coverage(ali_pafr, target=FALSE)
    expect_is(p_q, "ggplot")
    expect_equal(layer_data(p_q)[1,2],6273420)
})
