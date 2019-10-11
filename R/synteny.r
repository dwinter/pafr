#' Plot synteny between a query and target sequence in a paf alignment
#' 
#' @param ali pafr or tibble containing the genome alignment (as returned by
#' \code{\link{read_paf}})
#' @param q_chrom character, name for the query sequence
#' @param t_chrom character, name for the target sequence
#' @param rc logical, if TRUE, use the reverse and complement for the target
#' sequence
#' @param centre, logical, if TRUE (default) adjust the position of the target 
#' sequence so it is centred on the query. If not, both sequences start at
#' position zero
#' @param xlab, character, Name for for x-axis
#' @param ylab, character, Name for for x-axus
#' @param x_labeller function to be used to label x-axis (defaults to
#' @return A ggplot object dsplaying synteny between query and target sequences
#' @export
plot_synteny <- function(ali, q_chrom, t_chrom, centre = TRUE, rc = FALSE, 
                         xlab = "Position in query", ylab = "",
                         x_labeller=Mb_lab) {
    synt_df <- synteny_data(ali, q_chrom = q_chrom,  t_chrom = t_chrom, rc = rc)
    cs <- chrom_sizes(ali)
    tlen <- cs$tlens$tlen[cs$tlens$tname == t_chrom]
    qlen <- cs$qlens$qlen[cs$qlens$qname == q_chrom]
    seq_lens <- data.frame(seq = factor(c(q_chrom, t_chrom)),
                           start = 0, end = c(qlen, tlen))
    if (centre) {
        os <- abs(tlen - qlen) / 2
        if (qlen > tlen) {
            synt_df[["x"]][synt_df[["seq"]] == t_chrom ] <- synt_df[["x"]][synt_df[["seq"]]  == t_chrom ] + os
            seq_lens[2, 2:3] <- seq_lens[2, 2:3]  + os

        } else {
            synt_df[["x"]][synt_df[["seq"]] == q_chrom] <- synt_df[["x"]][ synt_df[["seq"]]  == q_chrom ] + os
            seq_lens[1, 2:3] <- seq_lens[1, 2:3]  + os
        }
    }
    if (q_chrom > t_chrom) {
        bottoms <- c(2.05, 0.8)
        tops <- c(2.2, 0.95)
    } else {
        bottoms <- c(0.8, 2.05)
        tops <- c(0.95, 2.2)
    }
    ggplot() + geom_polygon(data = synt_df, aes(x, seq, group = block_id),
                            fill = "grey80", colour = "black") +
               geom_rect(data = seq_lens,
                     aes(xmin = start, xmax = end, ymin = bottoms, ymax = tops),
                     colour = "black", fill = "white") +
               scale_x_continuous(xlab, labels = x_labeller) + ylab(ylab)
}

      
synteny_data <- function(ali, q_chrom, t_chrom, rc = FALSE) {
    to_plot <- subset(ali, qname == q_chrom & tname == t_chrom)
    if (rc) {
        e <- to_plot$tend
        s <- to_plot$tstart
        to_plot$tstart <- to_plot$tlen - e
        to_plot$tend <- to_plot$tlen - s
        to_plot$strand <- ifelse(to_plot$strand == "+", "-", "+")
    }
    synt_data <- apply(to_plot, 1, .per_row,
                       q_chrom = q_chrom, t_chrom = t_chrom)
    synt_df <- do.call(rbind.data.frame, synt_data)
    synt_df$block_id <- rep(seq_len(length(synt_data)), each = 5)
    synt_df
}

.per_row <- function(r, q_chrom, t_chrom) {
    if (r[["strand"]] == "+") {
        cols <- c("qstart", "qend", "tend", "tstart", "qstart")
    } else {
        cols <- c("qend", "qstart", "tend", "tstart", "qend")
    }
    list(x = as.numeric(r[cols]),
         seq = c(q_chrom, q_chrom, t_chrom, t_chrom, q_chrom))
}
