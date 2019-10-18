



#' Extract the sizes of all sequences in a paf alignment
#'
#' @param ali pafr or tibble containing the genome alignment (as returned by
#' \code{\link{read_paf}})
#' @return list with two elements (tlens and qlens). Each element is a
#' data.frame with one column of sequence names and another column containing
#' the length of each sequence
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package = "pafr") )
#' chrom_sizes(ali)
#' @export
chrom_sizes <- function(ali) {
    res <- list(qlens = unique(ali[, c("qname", "qlen")]),
                tlens = unique(ali[, c("tname", "tlen")]))
    if ("pafr" %in% class(ali)) {
        return(lapply(res, as.data.frame))
    }
    res
}


# Internal function used for concatenating sequences into a single ordering,
# used by dotplot and associated functions. Returns a list of four elements,
# which represent the starting positions of each sequence in a concatinated
# alignment and the total length of the query and targets sequences
#
#' @importFrom dplyr top_n
#' @importFrom dplyr group_by
#' @importFrom utils head
order_seqs <- function(ali, by, ordering = list()) {
    chrom_lens <- chrom_sizes(ali)
    qsum <- sum(chrom_lens[["qlens"]][, 2])
    tsum <- sum(chrom_lens[["tlens"]][, 2])
    if (by == "size") {
        q_idx <- order(chrom_lens[["qlens"]][, 2], decreasing = TRUE)
        t_idx <- order(chrom_lens[["tlens"]][, 2], decreasing = TRUE)
        qmap  <- structure(.Names = chrom_lens[["qlens"]][q_idx, 1],
                        c(0, head(cumsum(chrom_lens[["qlens"]][q_idx, 2]), -1)))
        tmap  <- structure(.Names = chrom_lens[["tlens"]][t_idx, 1],
                        c(0, head(cumsum(chrom_lens[["tlens"]][t_idx, 2]), -1)))
    } else if (by == "qstart") {
        #TODO
        #qidx/map id duplicated from above. DRY out depending on how we
        #implement other options
        q_idx <- order(chrom_lens[["qlens"]][, 2], decreasing = TRUE)
        qmap  <- structure(.Names = chrom_lens[["qlens"]][q_idx, 1],
                        c(0, head(cumsum(chrom_lens[["qlens"]][q_idx, 2]), -1)))
        longest_by_target <- top_n(group_by(ali, tname), 1, alen)
        t_idx <- order(qmap[longest_by_target$qname] + longest_by_target$qstart)
        tmap <- sort(structure(.Names = longest_by_target$tname[t_idx],
                         c(0, head(cumsum(longest_by_target$tlen[t_idx]), -1))))        
    } else if (by == "provided") {
        if (length(ordering) != 2) {
            stop("To use 'provided' sequence ordering argument 'ordering' must by a list with two character vectors")
        }
        qord <- ordering[[1]]
        tord <- ordering[[2]]
        q_idx <- match(qord, chrom_lens[["qlens"]][["qname"]])
        if (any(is.na(q_idx))) {
            msg <- paste("Sequence(s) provided for ordering are not present in alignment:\n", qord[is.na(q_idx)], "\n")
            stop(msg)
        }
        t_idx <- match(tord, chrom_lens[["tlens"]][["tname"]])
        if (any(is.na(t_idx))) {
            msg <- paste("Sequence(s) provided for ordering are not present in alignment:\n", tord[is.na(t_idx)], "\n")
            stop(msg)
        }
        qmap <- structure(.Names = chrom_lens[["qlens"]][q_idx, 1],
                        c(0, head(cumsum(chrom_lens[["qlens"]][q_idx, 2]), -1)))
        tmap <- structure(.Names = chrom_lens[["tlens"]][t_idx, 1],
                        c(0, head(cumsum(chrom_lens[["tlens"]][t_idx, 2]), -1)))
    }
    list(qmap = qmap, qsum = qsum, tmap = tmap, tsum = tsum)
}
               

add_pos_in_concatentaed_genome <- function(ali, maps) {
    ali$concat_qstart <- ali$qstart + maps[["qmap"]][ali$qname]
    ali$concat_qend <- ali$qend + maps[["qmap"]][ali$qname]
    ali$concat_tstart <- ali$tstart + maps[["tmap"]][ali$tname]
    ali$concat_tend <- ali$tend + maps[["tmap"]][ali$tname]
    ali
}


dotplot_name_df <- function(seq_map, genome_len) {
    data.frame(seq_name = names(seq_map),
               centre = seq_map + diff(c(seq_map, genome_len) / 2))
}


check_ordering <- function(ali, ordering) {
    q_in_order <- unique(ali[["qname"]]) %in% ordering[[1]]
    t_in_order <- unique(ali[["tname"]]) %in% ordering[[2]]
    if (any(!q_in_order)) {
        msg <- paste("Dropping data from sequences absent from ordering:\n",
                     paste(unique(ali[["qname"]])[!q_in_order], collapse = ","))
        warning(msg, call. = FALSE)
    }
    if (any(!t_in_order)) {
        msg <- paste("Dropping data from sequences absent from ordering:\n",
                     paste(unique(ali[["tname"]])[!t_in_order], collapse = ","))
        warning(msg, call. = FALSE)
    }
    return(invisible())
}

#' Generate a dotplot from a paf alignment.
#'
#' @param ali pafr or tibble containing the genome alignment (as returned by
#' \code{\link{read_paf}})
#' @param order_by How should the query and target sequences be ordered in the
#' dotplot. One of 'size' (smallest-to-largest), 'qstart' (query organised
#' smallest to largest, target by first match in the query genome or 'provided'
#' (ordering specified in \code{ordering} argument
#' @param label_seqs boolean, if \code{TRUE} label centre of query and target
#' seqeuences in margins of the dotplot
#' @param dashes boolen, if \code{TRUE} add dashes ad borders of query and
#' target sequences in the dotplot
#' @param ordering If \code{order_by} is set to \code{TRUE},
#' this should be set to a list with two elements, sepcifying the order of query
#' and then target sequences in the doptplot. This option is ignored if
#' \code{order_by} is set to other values
#' @param alignment_colour character, the colour used to draw each aligned
#' section in the dotplot (defaults to black).
#' @param xlab character The x-axis label (defaults to "query").
#' @param ylab character The y-axis label (defaults to "target").
#' @param line_size the width of the line used to represent an alignment in the
#' dotplot (defaults to 2)
#' @import ggplot2
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package = "pafr") )
#' dotplot(ali)
#' dotplot(ali) + theme_bw()
#' dotplot(ali, label_seqs=TRUE, order_by="qstart", alignment_colour="blue")
#' @export
dotplot <- function(ali, order_by = c("size", "qstart", "provided"),
                    label_seqs = FALSE, dashes = TRUE, ordering = list(),
                    alignment_colour = "black", xlab = "query", ylab = "target",
                    line_size=2) {
    by <- match.arg(order_by)
    if (by == "provided") {
        check_ordering(ali, ordering)
        ali <- subset(ali, qname %in% ordering[[1]] & tname %in% ordering[[2]])
    }
    seq_maps <- order_seqs(ali, by, ordering)
    ali <- add_pos_in_concatentaed_genome(ali, seq_maps)
    #now build the plot, first the aligments, easiest to do this +ve strand
    #first, then -ve
    p <- ggplot() +
      geom_segment(data = subset(ali, strand=="+"), 
         aes(x = concat_qstart, xend = concat_qend, y = concat_tstart, yend = concat_tend),
         size=line_size, colour=alignment_colour) +
      geom_segment(data = subset(ali, strand=="-"),
        aes(x = concat_qend, xend = concat_qstart, y = concat_tstart, yend = concat_tend),
        size = line_size, colour = alignment_colour) +
      coord_equal() +
      scale_x_continuous(xlab, labels = Mb_lab) + 
      scale_y_continuous(ylab, labels = Mb_lab)
  if (dashes) {
      p <- p + geom_hline(yintercept = c(seq_maps[["tmap"]], sum(unique(ali$tlen))), linetype = 3) +
               geom_vline(xintercept = c(seq_maps[["qmap"]], sum(unique(ali$qlen))), linetype = 3)
  }
  if (label_seqs) {
    qname_df <- dotplot_name_df(seq_maps[["qmap"]], seq_maps[["qsum"]])
    tname_df <- dotplot_name_df(seq_maps[["tmap"]], seq_maps[["tsum"]])
    p <- p + geom_text(data = qname_df, aes(label = seq_name, x = centre, y = 0),
                       vjust = 1, check_overlap = TRUE) +
             geom_text(data = tname_df,
                       aes(label = seq_name, x = 0, y = centre),
                       angle = 90, vjust = 0, check_overlap = TRUE)
  }
  p
}

highlight_dotplot <- function(hl_source, ali, bed, ordered_by, ordering, fill, colour, alpha){
    map_n <-  if (hl_source == "query") 1 else 3
    seq_maps <- order_seqs(ali, ordered_by)
    os <- seq_maps[[map_n]][as.character(bed[["chrom"]])]
    to_plot <- data.frame(i_start = bed$start + os, i_end = bed$end + os)
    if (hl_source == "query") {
        ystart <- 0
        yend <- seq_maps[["tsum"]]
        xstart <- "i_start"
        xend <- "i_end"
    } else {
        ystart <- "i_start"
        yend <- "i_end"        
        xend <- seq_maps[["qsum"]]
        xstart <- 0
    }

    geom_rect(data = to_plot,
             aes_string(xmin = xstart, xmax = xend, ymin = ystart, ymax = yend),
             fill = fill, colour = colour, alpha = alpha)
}

#' @rdname highlight_dotplot
#' @export 
highlight_query <- function(ali, bed, 
                           ordered_by = c("size", "qstart", "provided"), 
                           ordering = list(), fill = "yellow", colour = "black",
                           alpha=0.6) {
    highlight_dotplot("query", ali, bed,  
                      match.arg(ordered_by), ordering, fill, colour, alpha)
}
  
#' Highlight segments of a query or target genome in a dotplot
#'
#' This plot is intended to be used in conjunction with \code{link{dotplot}}.
#' Adding \code{higlight_query} or \code{highlight_target} to a dotplot call
#' (see examples below) will add a rectangular 'highlight' corresponding to a
#' particular genomic interval in the corresponding genome.
#'
#' @param ali pafr or tibble containing a genome alignment (as returned by
#' \code{\link{read_paf}})
#' @param bed \code{data.frame} or \code{tbl_df} containing a bed file, as returned by
#' \code{\link{read_bed}}. Should contain three columns, named "chrom", "start" 
#'  and "end".
#' @param ordered_by How are the query and target sequences be ordered in the
#' dotplot this highlight is being added to. This argument should match that
#' used in \code{\link{dotplot}}.
#' @param ordering If \code{order_by} is set to \code{TRUE},
#' this should be set to a list with two elements, sepcifying the order of query
#' and then target sequences in the doptplot. This option is ignored if
#' \code{order_by} is set to other values
#' @param fill character fill colour for highlight segment
#' @param colour, character outline colour for highlight segment
#' @param alpha, character opacity ([0-1]) for highlight segment
#' @rdname highlight_dotplot
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package = "pafr") )
#' cen <- read_bed(system.file("extdata", "Q_centro.bed", package = "pafr"))
#' dotplot(ali) + highlight_query(ali, cen)
#' interval <- data.frame(chrom="T_chr3", start = 2000000, end= 3000000)
#' dotplot(ali, label_seqs=TRUE) + 
#'    highlight_target(ali, interval)
#' @export
highlight_target <- function(ali, bed, 
                             ordered_by = c("size", "qstart", "provided"), 
                             ordering=list(), 
                             fill="yellow", colour="black", alpha=0.6){
    highlight_dotplot("target", ali, bed,  
                      match.arg(ordered_by), ordering, fill, colour, alpha)
}
