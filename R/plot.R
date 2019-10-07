
.per_row <- function(r, q_chrom, t_chrom){
    if(r[["strand"]] == "+"){
        cols <- c("qstart", "qend", "tend", "tstart", "qstart")
    } else {
        cols <- c("qend", "qstart", "tend", "tstart", "qend")
    }
    list(x =as.numeric( r[cols] ), seq = c(q_chrom, q_chrom, t_chrom, t_chrom, q_chrom))
}
      
synteny_data <- function(ali, q_chrom, t_chrom, RC=FALSE){
    to_plot <- subset(ali, qname==q_chrom & tname==t_chrom)
    if(RC){
        e <- to_plot$tend
        s <- to_plot$tstart
        to_plot$tstart = to_plot$tlen -e
        to_plot$tend = to_plot$tlen - s
        to_plot$strand <- ifelse(to_plot$strand == "+", "-", "+")
    }
    synt_data <- apply(to_plot, 1, .per_row, q_chrom=q_chrom, t_chrom=t_chrom)
    synt_df <- do.call(rbind.data.frame, synt_data)
    synt_df$block_id <- rep(1:length(synt_data), each=5)
    synt_df
}



highlight_dotpot <- function(hl_source, ali, bed, ordered_by, odering, fill, colour, alpha){
    stem <-  if (hl_source == "query") "q" else "t"
    seq_maps <- order_seqs(ali, ordered_by)
    os = seq_maps[[paste0(stem, "map")]][as.character(bed$chrom)]
    to_plot <- data.frame( i_start = bed$start + os, i_end = bed$end + os)
    if(hl_source == "query"){
        ystart = 0
        yend <- seq_maps[["tsum"]]
        xstart = "i_start"
        xend = "i_end"        
    } else {
        ystart = "i_start"
        yend = "i_end"        
        xend <- seq_maps[["qsum"]]
        xstart = 0
    }

    geom_rect(data=to_plot,
              aes_string(xmin=xstart, xmax=xend, ymin=ystart, ymax= yend),
               fill=fill, colour=colour, alpha=alpha)
}

#' @export 
highlight_query <- function(ali, bed, ordered_by = c("size", "qstart", "provided"), ordering=list(), fill="yellow", colour="black", alpha=0.6){
    by <- match.arg(ordered_by)
    seq_maps <- order_seqs(ali, by)
    os = seq_maps[["qmap"]][as.character(bed$chrom)]
    to_plot <- data.frame( qstart = bed$start + os, qend = bed$end + os)
    geom_rect(data=to_plot, 
              aes(xmin=qstart, xmax=qend, ymin=0, ymax=seq_maps[["tsum"]]),
              fill=fill, colour=colour, alpha=alpha)
}          
    
#' @export
highlight_target <- function(ali, bed, ordered_by = c("size", "qstart", "provided"), fill="yellow", colour="black", alpha=0.6){
    by <- match.arg(ordered_by)
    seq_maps <- order_seqs(ali, by)
    os = seq_maps[["tmap"]][as.character(bed$chrom)] 
    to_plot <- data.frame( tstart = bed$start + os, tend = bed$end + os) 
    geom_rect(data=to_plot, 
              aes(xmin=0, xmax=seq_maps[["qsum"]], ymin=tstart, ymax=tend), 
              fill=fill, colour=colour, alpha=alpha)
}



#' @export
chrom_sizes <- function(ali){
    res <- list(qlens = unique(ali[,c("qname", "qlen")]),
         tlens = unique(ali[,c("tname", "tlen")])
    )
    if('pafr' %in% class(ali)){
        return(lapply(res, as.data.frame))
    }
    res
}


#' Plot synteny between a query and target sequence in a paf alignment
#' 
#' @param ali pafr or tibble containing the genome alignment (as returned by
#' \code{\link{read_paf}})
#' @param q_chrom character, name for the query sequence
#' @param t_chrom character, name for the target sequence
#' @param RC logical, if TRUE, use the reverse and complement for the target
#' sequence
#' @param centre, logical, if TRUE (default) adjust the position of the target 
#' sequence so it is centred on the query. If not, both sequences start at
#' position zero
#' @param xlab, character, Name for for x-axis
#' @param ylab, character, Name for for x-axus
#' @param x_labeller function to be used to label x-axis (defaults to
#' @return A ggplot object dsplaying synteny between query and target sequences
#' @export
plot_synteny <- function(ali, q_chrom, t_chrom, centre=TRUE, RC=FALSE, 
                         xlab = "Position in query", ylab = "",  x_labeller=Mb_lab){
    synt_df <- synteny_data(ali, q_chrom = q_chrom,  t_chrom = t_chrom, RC=RC)
    cs <- chrom_sizes(ali)
    tlen <- cs$tlens$tlen[cs$tlens$tname == t_chrom]     
    qlen <- cs$qlens$qlen[cs$qlens$qname == q_chrom]
    seq_lens <- data.frame(seq=factor(c(q_chrom, t_chrom)), start=0, end=c(qlen, tlen))
    if( centre ){
        os <- abs(tlen - qlen)/2
        if(qlen > tlen){        
            synt_df[["x"]] [ synt_df[["seq"]]  == t_chrom ] <- synt_df[["x"]] [ synt_df[["seq"]]  == t_chrom ] + os
            seq_lens[2,2:3] <- seq_lens[2,2:3]  + os

        } else {
            synt_df[["x"]] [ synt_df[["seq"]]  == q_chrom ] <- synt_df[["x"]] [ synt_df[["seq"]]  == q_chrom ] + os
            seq_lens[1,2:3] <- seq_lens[1,2:3]  + os
        }
    }
    if(q_chrom > t_chrom){
        bottoms <- c(2.05, 0.8)
        tops <- c(2.2, 0.95)
    } else {
        bottoms <- c(0.8, 2.05)
        tops <- c(0.95, 2.2)
    }
    ggplot() + geom_polygon(data=synt_df, aes(x,seq, group=block_id), fill='grey80', colour='black') + 
               geom_rect(data=seq_lens, aes(xmin=start, xmax=end, ymin=bottoms, ymax=tops), colour='black', fill='white') +
               scale_x_continuous(xlab, label=x_labeller) + ylab(ylab)


} 


#' @importFrom dplyr top_n
#' @importFrom dplyr group_by
#' @importFrom utils head
order_seqs <- function(ali, by, ordering = list() ){      
    chrom_lens <- chrom_sizes(ali)
    qsum = sum(chrom_lens[["qlens"]][,2])
    tsum = sum(chrom_lens[["tlens"]][,2])
    if(by == "asis"){
        ("NOT YET IMPLEMENTED")
    }
    if(by == "size"){
        q_idx <- order(chrom_lens[["qlens"]][,2], decreasing = TRUE )
        t_idx <- order(chrom_lens[["tlens"]][,2], decreasing = TRUE )
        qmap  = structure(.Names=chrom_lens[["qlens"]][q_idx,1], c(0, head(cumsum(chrom_lens[["qlens"]][q_idx,2]),-1)))
        tmap  = structure(.Names=chrom_lens[["tlens"]][t_idx,1], c(0, head(cumsum(chrom_lens[["tlens"]][t_idx,2]),-1)))
    } else if (by == "qstart") {
        #TODO
        #qidx/map id duplicated from above. DRY out depending on how we
        #implement other options
        q_idx <- order(chrom_lens[["qlens"]][,2], decreasing = TRUE )
        qmap  = structure(.Names=chrom_lens[["qlens"]][q_idx,1], c(0, head(cumsum(chrom_lens[["qlens"]][q_idx,2]),-1)))
        longest_by_target <- top_n(group_by(ali, tname), 1, alen)
        t_idx  <- order(qmap[longest_by_target$qname] + longest_by_target$qstart)
        tmap <- sort(structure(.Names = longest_by_target$tname[t_idx], c(0, head(cumsum(longest_by_target$tlen[t_idx]), -1))))        
    } else if (by == "provided") {
        if(length(ordering) != 2) {
            stop("To use 'provided' sequence ordering argument 'ordering' must by a list with two character vectors")             
        }
        qord <- ordering[[1]]
        tord <- ordering[[2]]
        q_idx <- match(qord, chrom_lens[["qlens"]][["qname"]])
        if(any(is.na(q_idx))){
            msg <- paste("Sequence(s) provided for ordering are nor present in alignment:\n",qord[is.na(q_idx)], "\n")
            stop(msg)
        }
        t_idx <- match(tord, chrom_lens[["tlens"]][["tname"]])
        if(any(is.na(t_idx))){
            msg <- paste("Sequence(s) provided for ordering are nor present in alignment:\n",tord[is.na(t_idx)], "\n")
            stop(msg)
        }
        qmap  = structure(.Names=chrom_lens[["qlens"]][q_idx,1], c(0, head(cumsum(chrom_lens[["qlens"]][q_idx,2]),-1)))
        tmap  = structure(.Names=chrom_lens[["tlens"]][t_idx,1], c(0, head(cumsum(chrom_lens[["tlens"]][t_idx,2]),-1)))
    }
    list(qmap = qmap, qsum = qsum, tmap = tmap, tsum = tsum)
}
               

add_pos_in_concatentaed_genome <- function(ali, maps){
    ali$concat_qstart <- ali$qstart + maps[["qmap"]][ali$qname]
    ali$concat_qend <- ali$qend + maps[["qmap"]][ali$qname]
    ali$concat_tstart <- ali$tstart + maps[["tmap"]][ali$tname]
    ali$concat_tend <- ali$tend + maps[["tmap"]][ali$tname]
    ali
}


dotplot_name_df <- function(seq_map, genome_len){
    data.frame(seq_name = names(seq_map),
               centre=seq_map + diff(c(seq_map, genome_len)/2)
    )
}

check_ordering <- function(ali, ordering){
    q_in_order <- unique(ali[["qname"]]) %in% ordering[[1]]
    t_in_order <- unique(ali[["tname"]]) %in% ordering[[2]]
    if(any(!q_in_order)){
        msg <- paste("Dropping data from sequences absent from ordering:\n",
                     paste(unique(ali[["qname"]])[!q_in_order], collapse=","))
        warning(msg, call.=FALSE)
    }
    if(any(!t_in_order)){
        msg <- paste("Dropping data from sequences absent from ordering:\n",
                     paste(unique(ali[["tname"]])[!t_in_order], collapse=","))
        warning(msg, call.=FALSE)
    }
    return(invisible()) 
}

#' @import ggplot2 
#' @export
dotplot <- function(ali, order_by = c("size", "qstart", "provided", "asis"), label_seqs = FALSE, dashes=TRUE, ordering = list(), alignment_colour="black", xlab = "query", ylab="target", line_size=2){
    by <- match.arg(order_by)
    seq_maps <- order_seqs(ali, by, ordering)
    if(by == "provided"){
        check_ordering(ali, ordering)
        ali <- subset(ali, qname %in% ordering[[1]] & tname %in% ordering[[2]])
    }
    ali <- add_pos_in_concatentaed_genome(ali, seq_maps)
  
    
  p <- ggplot() + 
      geom_segment(data=subset(ali, strand=="+"), aes(x=concat_qstart, xend=concat_qend, y=concat_tstart, yend=concat_tend), size=line_size, colour=alignment_colour) + 
      geom_segment(data=subset(ali, strand=="-"), aes(x=concat_qend, xend=concat_qstart, y=concat_tstart, yend=concat_tend), size=line_size, colour=alignment_colour) + 
      coord_equal() +
      scale_x_continuous(xlab, labels=Mb_lab) + 
      scale_y_continuous(ylab, labels=Mb_lab)
  if(dashes){
      p <- p +  geom_hline(yintercept=c(seq_maps[["tmap"]], sum(unique(ali$tlen))), linetype=3) +
                geom_vline(xintercept=c(seq_maps[["qmap"]], sum(unique(ali$qlen))), linetype=3)
  }
  if (label_seqs){
    qname_df <- dotplot_name_df(seq_maps[["qmap"]], seq_maps[["qsum"]])
    tname_df <- dotplot_name_df(seq_maps[["tmap"]], seq_maps[["tsum"]])
    p <- p + geom_text(data=qname_df, aes(label=seq_name, x=centre, y=0), vjust=1, check_overlap=TRUE) + 
             geom_text(data=tname_df, aes(label=seq_name, x=0, y=centre), angle=90, vjust=0, check_overlap=TRUE) 
  }
  p

}
