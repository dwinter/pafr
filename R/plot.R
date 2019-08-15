
.per_row <- function(r, q_chrom, t_chrom){
    cols 
    if(r[["strand"]] == "+"){
        cols <- c("qstart", "qend", "tend", "tstart", "qstart")
    } else {
        cols <- c("qend", "qstart", "tend", "tstart", "qend")
    }
    list(x =as.numeric( r[cols] ), seq = c(q_chrom, q_chrom, t_chrom, t_chrom, q_chrom))
}
      
#' @export
synteny_data <- function(ali, q_chrom, t_chrom, offset=0){
    to_plot <- subset(ali, qname==q_chrom & tname==t_chrom)
    synt_data <- apply(to_plot, 1, .per_row, q_chrom=q_chrom, t_chrom=t_chrom)
    synt_df <- do.call(rbind.data.frame, synt_data)
    synt_df$block_id <- rep(1:length(synt_data), each=5)
    synt_df
}


#' @export
Mb_lab <- function(x) paste(x/1e6, "Mb")

chrom_sizes <- function(ali){
    list(qlens = unique(ali[,c("qname", "qlen")]),
         tlens = unique(ali[,c("tname", "tlen")])
    )
}

#' @export
plot_synteny <- function(ali, q_chrom, t_chrom, centre=TRUE){
    synt_df <- synteny_data(ali, q_chrom = q_chrom,  t_chrom = t_chrom)
    cs <- chrom_sizes(ali)
    tlen <- cs$tlens$tlen[cs$tlens$tname == t_chrom]     
    qlen <- cs$qlens$qlen[cs$qlens$qname == q_chrom]
    os <- abs(tlen - qlen)/2
    seq_lens <- data.frame(seq=factor(c(q_chrom, t_chrom)), start=0, end=c(qlen, tlen))
    if( centre ){
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
               geom_rect(data=seq_lens, aes(xmin=start, xmax=end, ymin=bottoms, ymax=tops), colour='black', fill='white')

} 

#' @export
plot_coverage <- function(ali, target=TRUE, fill_colour="forestgreen"){
    cs <- chrom_sizes(ali)
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
    biggest <- max(u_chroms$seq_len)
    ggplot() + 
        geom_rect(data=u_chroms, aes(xmin=0, xmax=seq_len, ymin=-1, ymax=1), colour='black', fill='white', size=1) + 
        geom_rect(data=ali, aes_string(xmin=spos, xmax=epos, ymin=-0.95, ymax=0.95), fill=fill_colour) + 
        geom_text(data=u_chroms, aes(x=biggest/10,y=0, label=seq_name)) + 
        facet_grid(seq_name ~ . )  +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) 
    
}
    


order_seqs <- function(ali){
    t_idx <- order(unique(ali$tlen), decreasing=TRUE)
    q_idx <- order(unique(ali$qlen), decreasing=TRUE)
    list(qmap  = structure(.Names=unique(ali$qname)[q_idx], c(0, head(cumsum(unique(ali$qlen)[q_idx]),-1))),
         tmap  = structure(.Names=unique(ali$tname)[t_idx], c(0, head(cumsum(unique(ali$tlen)[t_idx]),-1)))
    )
               
}

add_pos_in_concatentaed_genome <- function(ali, maps){
    ali$concat_qstart <- ali$qstart + maps[["qmap"]][ali$qname]
    ali$concat_qend <- ali$qend + maps[["qmap"]][ali$qname]
    ali$concat_tstart <- ali$tstart + maps[["tmap"]][ali$tname]
    ali$concat_tend <- ali$tend + maps[["tmap"]][ali$tname]
    ali
}



#' @import ggplot2 
#' @export
dotplot <- function(ali, dashes=TRUE, alignment_colour="black", xlab = "query", ylab="target"){
  seq_maps <- order_seqs(ali) 
  ali <- add_pos_in_concatentaed_genome(ali, seq_maps)
  
    
  p <- ggplot() + 
      geom_segment(data=subset(ali, strand=="+"), aes(x=concat_qstart, xend=concat_qend, y=concat_tstart, yend=concat_tend), size=2, colour=alignment_colour) + 
      geom_segment(data=subset(ali, strand=="-"), aes(x=concat_qend, xend=concat_qstart, y=concat_tstart, yend=concat_tend), size=2, colour=alignment_colour) + 
      coord_equal() +
      scale_x_continuous(xlab, labels=Mb_lab) + 
      scale_y_continuous(ylab, labels=Mb_lab)
  if(dashes){
      return(p + geom_hline(yintercept=seq_maps[["tmap"]], linetype=3) +
                 geom_vline(xintercept=seq_maps[["qmap"]], linetype=3))
  }
  p

}
