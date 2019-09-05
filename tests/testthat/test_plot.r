context("Testing plotting... to the extent possible")

ali_pafr <- read_paf("test_ali.paf")

test_that("dotplot works produces a plot", {
    p <- dotplot(ali_pafr)
    expect_is(p, "ggplot")
    labs <- unlist(p[["labels"]])
    expect_equal(names(labs), c("x", "y", "xend", "yend", "yintercept", "xintercept"))
    expect_equal(unname(labs), c("concat_qstart", "concat_tstart", "concat_qend", "concat_tend", "yintercept", "xintercept"))
})

test_that("Ordering sequences by size", {
    seq_map <- order_seqs(ali_pafr)
    qmap <- seq_map[["qmap"]]
    tmap <- seq_map[["tmap"]]
    expect_equal(names(qmap), c("Q_chr1", "Q_chr2", "Q_chr6", "Q_chr7"))
    expect_equal(unname(qmap), c(0, 6273420, 12125925, 16469409))
    expect_equal(names(tmap), c("T_chr1", "T_chr2", "T_chr3", "T_chr5"))
    expect_equal(unname(tmap), c(0, 5893594, 10514309, 16054194))
})

test_that("Ordering sequences by size", {
    seq_map <- order_seqs(ali_pafr, "qstart")
    tmap <- seq_map[["tmap"]]
    expect_equal(names(tmap), c("T_chr2", "T_chr1", "T_chr5", "T_chr3"))
    expect_equal(unname(tmap), c(0, 5893594, 12271384, 16892099))
})

test_that("Ordering sequences to preference", {
    chrom_order <-  list(
        c("Q_chr6", "Q_chr1", "Q_chr2", "Q_chr2", "Q_chr7"), 
        c("T_chr5", "T_chr3", "T_chr2", "T_chr3", "T_chr1")
    )              
    seq_map <- order_seqs(ali_pafr, "provided", chrom_order)
    qmap <- seq_map[["qmap"]]
    tmap <- seq_map[["tmap"]]
    expect_equal(names(qmap), chrom_order[[1]])
    expect_equal(names(tmap), chrom_order[[2]])
})

test_that("Ordering sequences to preference errors out with no chroms", {
    expect_error(seq_map <- order_seqs(ali_pafr, "provided"))
}) 

test_that("Ordering sequences to preference erros with missing chroms", {
    chrom_order <-  list(
        c("Q_chr6", "Q_chr1", "Q_chr2", "Q_chr2", "Q_chr7", "another one"), 
        c("T_chr5", "T_chr3", "T_chr2", "T_chr3", "T_chr1")
    )              
    expect_error(seq_map <- order_seqs(ali_pafr, "provided", chrom_order))
})

test_that("dotplot options change plot layers", {
    p <- dotplot(ali_pafr)
    p_lab <- dotplot(ali_pafr, label_seqs=TRUE)
    p_nodash <- dotplot(ali_pafr, dashes=FALSE)
    expect_length(p[["layers"]], 4)
    expect_length(p_lab[["layers"]], 6)
    expect_length(p_nodash[["layers"]], 2)
})
