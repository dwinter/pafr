#' Plot the regions of one genome that are covered by alignments in a paf file
#'
#' Each sequence in the focal genome is displayed as a rectangle, with regions
#' covered by an alignment shaded as per the \code{fill} argument described
#' below. Uncovered regions remain white.
#'
#' Note that this function uses \code{\link{theme_coverage_plot}} to style the
#' graph. Using another ggplot theme on the plot may produce unexpected results.
#'
#' @param ali alignment  An alignment as read by \code{\link{read_paf}}
#' @param target logical  If TRUE, display coverage for the target
#' genome; if FALSE, display coverage for the query
#' @param fill character  How to colour the alignment blocks. If the character
#' provided is the name of a column in the alignment, this column will be passed
#' to \code{\link{ggplot2}} to shade alignment blocks. Otherwise, the character
#' is treated as a single colour to be used for all alignment blocks.
#' @param direct_label logical  If TRUE, use geom_text to directly label the
#' name of the focal sequences; if FALSE, no direct labels are drawn
#' @param label_colour character  Colour used for direct labels
#' @param xlab string  Name for the x-axis
#' @param x_labeller function  Function to be used to label the x-axis (defaults to
#' \code{\link{Mb_lab}}
#' @importFrom rlang .data
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package="pafr") )
#' plot_coverage(ali)
#' plot_coverage(ali, fill='qname', direct_label=FALSE) + 
#'    scale_fill_brewer(palette="Set1")
#' @export
plot_coverage <- function(ali, target = TRUE, fill = "forestgreen",
                          direct_label = TRUE, label_colour = "black",
                          xlab = "Position in sequence", x_labeller = Mb_lab) {
    cs <- chrom_sizes(ali)
    #first ID the appropriate columsns to use in the plot (target or query
    #seqeunces and the locatoins of alignments on those sequences).
    if (target) {
        u_chroms <- cs[["tlens"]]
        spos <- "tstart"
        epos <- "tend"
        names(ali)[6] <- "seq_name"
    } else {
        u_chroms <- cs[["qlens"]]
        spos <- "qstart"
        epos <- "qend"
        names(ali)[1] <- "seq_name"
    }
    names(u_chroms) <- c("seq_name", "seq_len")
    #start building the plot, always contains sequence outlines:

    p <- ggplot() +
        geom_rect(data = u_chroms,
                  aes(xmin = 0, xmax = seq_len, ymin = -1, ymax = 1),
                  colour = "black", fill = "white", size = 1)

    if (fill %in% colnames(ali)) {
        #use column as an attribute for plot
        p <- p + geom_rect(data = ali, aes_string(xmin = spos, xmax = epos, ymin = -0.95, ymax = 0.95, fill = fill))
    } else {
        #If fill is not a column name, assume it is a single colour to use for
        #all algignments
        p <- p + geom_rect(data = ali, aes_string(xmin = spos, xmax = epos, ymin = -0.95, ymax = 0.95), fill = fill)
    }
    #However we make the plot, it needs some polish:
    p <- p + scale_x_continuous(xlab, labels = x_labeller) +
        facet_grid(.data[["seq_name"]] ~ .)
    #finally label it one way or another and return.
    if (direct_label) {
        biggest <- max(u_chroms$seq_len)
        p <- p + geom_text(data = u_chroms,
                           aes(x = biggest / 10,y = 0, label = .data[["seq_name"]]),
                           color = label_colour) +
            theme_coverage_plot(facet_labs=FALSE)
        return(p)
    }
    p  + theme_coverage_plot()
}

#' A minimalistic ggplot2 theme designed for use with genome coverage plots
#'
#' This theme is used as the default when \code{\link{plot_coverage}} is called,
#' so you should usually only call this function to modify the appearance of the
#' coverage plot.
#' @param facet_labs logical  If TRUE (default), label sequences using the facet
#' labels; if FALSE, sequences are labeled directly using 
#' \code{\link[ggplot2]{geom_text}}
#' @param show_legend logical  If TRUE (default), label display any legend
#' associated with the fill parameter of \code{plot_coverage}; if FALSE, 
#' do not display a legend
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package="pafr") )
#' plot_coverage(ali) + theme_coverage_plot(show_legend=FALSE)
#' @export
theme_coverage_plot <- function(facet_labs = TRUE, show_legend = TRUE) {
    theme <- theme_bw()
    theme$axis.title.y <- element_blank()
    theme$axis.text.y <- element_blank()
    theme$axis.ticks.y <- element_blank()
    if (!facet_labs) {
          theme$strip.background <- element_blank()
          theme$strip.text <- element_blank()
    }
    if (!show_legend) {
        theme$legend.position <- "none"
    }
    theme
}
