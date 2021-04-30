write_simple_table <- function (...) {
    write.table(..., quote = FALSE, row.names = FALSE, col.names = FALSE, 
        sep = "\t")
}



circos_ktype <- function(paf_table){
    cs <- chrom_sizes(paf_table)
    q_sizes <- cs[["qlens"]]
    t_sizes <- cs[["tlens"]]
    q_ktype <- data.frame(type = "chr", parent_structure = "-",                           
                           ID = q_sizes[["qname"]], label = 1:nrow(q_sizes),
                           start = 0, end = q_sizes[["qlen"]],      
                           colour = paste("chr", 1:nrow(q_sizes), sep = "")
    )
    t_ktype <- data.frame(type = "chr", parent_structure = "-",
                           ID = t_sizes[["tname"]], label = 1:nrow(t_sizes),
                           start = 0, end = t_sizes[["tlen"]],      
                           colour = paste("chr", 1:nrow(t_sizes), sep = "")
    )
    return(list (query = q_ktype, target = t_ktype))
}

circos_links <- function(paf_table){
    res <- paf_table[,c(1,3,4,6,8,9)]
    class(res) <- c("data.frame")
    res

}

paf_to_circos <- function(paf_table, file_stem=NULL, return =TRUE){
    ktypes <- circos_ktype(paf_table)
    links <- circos_links(paf_table)
    if(!is.null(file_stem)){
        write_simple_table(ktypes[["query"]], paste0(file_stem, "target.karyotype"))
        write_simple_table(ktypes[["target"]], paste0(file_stem, "query.karyotype"))
        write_simple_table(links, paste0(file_stem, "links.txt"))
    }
    if(return){
        return( c(ktypes, links=links))
    }

}

#usage
#
# fungi <- read_paf("inst/extdata/fungi.paf")
# clean_fungi <- filter_secondary_alignments(fungi)
# clean_fungi <- subset(clean_fungi, alen > 1e4 & mapq > 40)
# write.simple.table
# paf_to_circos(clean_fungi, file_stem="fungi_", return=FALSE)
#
