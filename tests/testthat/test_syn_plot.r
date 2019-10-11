
context("Test the sytnteny plots work")

test_that("plot_syntent produces a plot", {
    p <- plot_synteny(ali_pafr, "Q_chr1", "T_chr2")
    expect_is(p, "ggplot")

})

test_that("synteny data can be extracted", {
    d1 <- synteny_data(ali_pafr, "Q_chr1", "T_chr2")
    expect_equal(d1[["x"]], c(665479, 660794, 5228905, 5224220, 665479))
})

test_that("RC flips only the query chrom", {
    d2 <- synteny_data(ali_pafr, "Q_chr1", "T_chr2", rc = TRUE)
    expect_equal(d2[["x"]], c(660794, 665479, 669374, 664689, 660794))
})


test_that("Centering moves the target chrom", {
    p <- plot_synteny(ali_pafr, "Q_chr1", "T_chr2")
    p2 <- plot_synteny(ali_pafr, "Q_chr1", "T_chr2", centre=FALSE)
    old_x <-  p[["layers"]][[1]][["layer_data"]]()[["x"]]
    new_x <- p2[["layers"]][[1]][["layer_data"]]()[["x"]]
    expect_equal(new_x, c(665479, 660794, 5228905, 5224220, 665479))
    expect_equal(old_x, c(665479, 660794, 5418818, 5414133, 665479))
})



