#' Plot synteny between a query and target sequence in a PAF alignment
#' 
#' @param ali pafr or tibble containing the genome alignment (as returned by
#' \code{\link{read_paf}})
#' @param q_chrom character  Name for the query sequence
#' @param t_chrom character  Name for the target sequence
#' @param rc logical  If TRUE, use the reverse and complement for the target
#' sequence
#' @param centre logical  If TRUE (default), adjust the position of the target 
#' sequence, so it is centred on the query. If not, both sequences start at
#' position zero
#' @param xlab string  Name for the x-axis
#' @param ylab string  Name for the y-axis
#' @param x_labeller  Function to be used to label the x-axis
#' @return  A ggplot object that displays synteny between query and target sequences
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package="pafr") )
#' long_ali <- subset(ali, alen > 1e4)
#' plot_synteny(long_ali, q_chrom="Q_chr3", t_chrom="T_chr4", centre=TRUE)
#' plot_synteny(long_ali, q_chrom="Q_chr5", t_chrom="T_chr5", centre=TRUE)
#' plot_synteny(long_ali, q_chrom="Q_chr5", t_chrom="T_chr5", centre=TRUE, rc=TRUE)
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
    ggplot() + geom_polygon(data = synt_df, aes_string(x="x", y="seq", group = "block_id"),
                            fill = "grey80", colour = "black") +
               geom_rect(data = seq_lens,
                     aes_string(xmin = "start", xmax = "end", ymin = "bottoms", ymax = "tops"),
                     colour = "black", fill = "white") +
               scale_x_continuous(xlab, labels = x_labeller) + ylab(ylab)
}

# internal function to convert alignments into polygon co-ordinates for plotting      
synteny_data <- function(ali, q_chrom, t_chrom, rc = FALSE) {
    to_plot <- ali[ali$qname == q_chrom & ali$tname == t_chrom,]
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
