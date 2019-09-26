
#' Plot the regions of one genome that are covered by alignments in a paf flie
#' 
#' Each sequence in the focal genome is display as rectange, with regions 
#' covered by an alignment shaded with the fill_colour as described below.
#' Uncovered regions remain white.
#'
#' Note this function uses a \code{\link{theme_coverage_plot}} to style the graph
#' using another ggplot theme on the plot may produce unexpected results.
#'
#' @param ali an alignment as read by \code{\link{read_paf}}
#' @param target logical, if TRUE (defautl) dsiplay coverage for the target
#' genome, if FALSE for the query
#' @param fill_colour, colour used for the covered-regions
#' @param direct_label, logical. If TRUE use geom_text to directly label the
#' name of the focal sequences. if FALSE no direct labels are drown
#' @param label_colour character, colour used for direct labels
#' @param xlab, character, Name for for x-axis
#' @param x_labeller function to be used to label x-axis (defaults to
#' \code{\link{Mb_lab}}
#' @export
plot_coverage <- function(ali, target=TRUE, fill="forestgreen", 
                          direct_label = TRUE, label_colour="black",
                          xlab="Position in sequence", x_labeller = Mb_lab){


    cs <- chrom_sizes(ali)
    #first ID the appropriate columsns to use in the plot (target or query
    #seqeunces and the locatoins of alignments on those sequences).    
    if(target){
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
        geom_rect(data=u_chroms, aes(xmin=0, xmax=seq_len, ymin=-1, ymax=1), colour='black', fill='white', size=1) 
    
    if(fill %in% colnames(ali)){
        #use column as an attribute for plot
        p <- p + geom_rect(data=ali, aes_string(xmin=spos, xmax=epos, ymin=-0.95, ymax=0.95, fill=fill)) 
    } else {
        #If fill is not a column name, assume it is a single colour to use for
        #all algignments
        p <- p + geom_rect(data=ali, aes_string(xmin=spos, xmax=epos, ymin=-0.95, ymax=0.95), fill=fill) 
    }
    #However we make the plot, it needs some polish:
    p <- p + scale_x_continuous(xlab, labels=x_labeller) +
        facet_grid(seq_name ~ . )  
    #finally label it one way or another and return.   
    if (direct_label){
        biggest <- max(u_chroms$seq_len)
        p <- p + geom_text(data=u_chroms, aes(x=biggest/10,y=0, label=seq_name), color=label_colour) +
            theme_coverage_plot(facet_labs=FALSE)
        return(p)
    }    
    p  + theme_coverage_plot() 
}
    
#' A minimalistic ggplot2 theme designed for use with genome coverage plots
#' @param facet_labs, logical. If TRUE (default) label sequences using the facet
#' labes. If FALSE sequences are labeled directly using \code{\link{geom_text}}
#' @export
theme_coverage_plot <- function(facet_labs=TRUE){
    theme <- theme_bw()
    theme$axis.title.y=element_blank()
    theme$axis.text.y=element_blank()
    theme$axis.ticks.y=element_blank()
    if(!facet_labs){
          theme$strip.background = element_blank()
          theme$strip.text = element_blank()
    }
    theme
}


