context("Testing plotting... to the extent possible")

test_that("dotplot works produces a plot", {
    p <- dotplot(ali_pafr)
    expect_is(p, "ggplot")
    labs <- unlist(p[["labels"]])
    expect_equal(unname(labs["x"]), "concat_qstart")
    expect_equal(unname(labs["xend"]), "concat_qend")
    expect_equal(unname(labs["y"]), "concat_tstart")
    expect_equal(unname(labs["yend"]), "concat_tend")
    expect_equal(unname(labs["xintercept"]), "xintercept")
    expect_equal(unname(labs["yintercept"]), "yintercept")
})

test_that("Ordering sequences by size", {
    seq_map <- order_seqs(ali_pafr, by="size")
    qmap <- seq_map[["qmap"]]
    tmap <- seq_map[["tmap"]]
    expect_equal(names(qmap), c("Q_chr1", "Q_chr2", "Q_chr6", "Q_chr7"))
    expect_equal(unname(qmap), c(0, 6273420, 12125925, 16469409))
    expect_equal(names(tmap), c("T_chr1", "T_chr2", "T_chr3", "T_chr5"))
    expect_equal(unname(tmap), c(0, 6377790, 12271384, 17811269))
})

test_that("Ordering sequences by start", {
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

test_that("Ordering sequences to preference warns on drop", {
    chrom_order <-  list(
        c("Q_chr6", "Q_chr1"), 
        c("T_chr3", "T_chr2")
    )
    expect_warning(dotplot(ali_pafr, "provided", ordering=chrom_order))
})

test_that("dotplot options change plot layers", {
    p <- dotplot(ali_pafr)
    p_lab <- dotplot(ali_pafr, label_seqs=TRUE)
    p_nodash <- dotplot(ali_pafr, dashes=FALSE)
    expect_length(p[["layers"]], 4)
    expect_length(p_lab[["layers"]], 6)
    expect_length(p_nodash[["layers"]], 2)
})

test_that("Can add query highlights to dotplots", {
    B <- read_bed("test_I.bed")
    p <- dotplot(ali_pafr) + highlight_query(B)
    #previously tested defaul has four layers
    expect_length(p[["layers"]], 5)
    mapping <- sapply(p[["layers"]][[5]][["mapping"]], quo_name)
    expect_equal(unname(mapping), c("istart", "iend", "0", "len"))
    expect_equal(names(mapping), c("xmin", "xmax", "ymin", "ymax"))
    
})

test_that("Can add target highlights to dotplots", {
    B <- read_bed("test_T.bed")
    p <- dotplot(ali_pafr) + highlight_target(B)
    expect_length(p[["layers"]], 5)
    mapping <- sapply(p[["layers"]][[5]][["mapping"]], quo_name)
    expect_equal(unname(mapping), c("0", "len", "istart", "iend"))
    expect_equal(names(mapping), c("xmin", "xmax", "ymin", "ymax"))
}) 

test_that("dotplot highlight warnings/errors", {
    B <- read_bed("test_T.bed")
    expect_error( dotplot(ali_pafr) + highlight_query(B))
    B$chrom[1] <- "X"
    expect_warning( dotplot(ali_pafr) + highlight_target(B))

})
