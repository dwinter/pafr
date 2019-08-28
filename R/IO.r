.make_numeric <- function(df, cols){
    for( i in cols){
        df[,i] <- as.numeric(df[,i])
    }
    df
}

tag_to_type_fxn <- function(tag_tokens){
    if(tag_tokens[[2]] %in% c("f", "H", "i")){
        return(as.numeric)
    }
    identity
}      


.one_tag_row <- function(tags){
    tokens <- strsplit(tags, ":")
    F <- sapply(tokens , tag_to_type_fxn)
    structure(lapply(1:(length(tokens)), 
           function(i) F[[i]](tokens[[i]][3])), .Names=sapply(tokens, "[[", 1))
    
}

#' @importFrom dplyr bind_rows
process_tags <- function(tag_rows){
    tag_tibble <- dplyr::bind_rows(lapply(tag_rows, .one_tag_row))
    as.data.frame(tag_tibble) #base format is df, so force it back here
}


#' @importFrom tibble as_tibble 
#' @export
read_paf <- function(file_name, tibble=FALSE){
    lines <- scan(file_name, "", sep="\n", quiet=TRUE)
    tokens <-  strsplit(lines, "\t")
    res <- do.call(rbind.data.frame, c(lapply(tokens, "[", 1:12), stringsAsFactors=FALSE))
    names(res) <- c("qname", "qlen", "qstart", "qend", "strand", 
                          "tname", "tlen", "tstart", "tend", "nmatch", 
                          "alen", "mapq")

    res <- .make_numeric(res, c(2,3,4,7:12))
    if(any(lengths(tokens) > 12)){
        raw_tags <- sapply(tokens, function(x) paste(x[13:length(x)]))
        res <- cbind.data.frame(res, process_tags(raw_tags))
    }
    if(tibble){
        return(as_tibble(res))
    }
    class(res) <- c("pafr", "data.frame")
    res
}

#' @export 
print.pafr <- function(x, ...){
    if(nrow(x) == 1){
        print_alignment(x)
        return(invisible(x))
    } 
    print_pafr_rows(x)
    return(invisible(x))
}
    
print_pafr_rows <- function(x){
    cat("pafr object with ", nrow(x), " alignments (", round(sum(x$alen)/1e6,1), "Mb)\n", sep="")
    cat(" ", length(unique(x$qname)), " query seqs\n", sep="")
    cat(" ", length(unique(x$tname)), " target seqs\n", sep="")
    ntags <- ncol(x) - 12
    if(ntags < 12){
        if(ntags > 0){
            tags <- paste(names(x)[13:ncol(x)], collapse=", ")
        } else {
            tags <- ""
        }

    } else {
        tags <- paste(paste(names(x)[13:23], collapse=", "), "...")
    }
    cat(" ", ntags, " tags: ", tags, "\n", sep="")
}

print_alignment <- function(x){
    ali_str <- paste0(x[["qname"]], ":", x[["qstart"]], "-", x[["qend"]], " v ", x[["tname"]], ":", x[["tstart"]], "-", x[["tend"]])
    cat("Single pafr alignment:\n")
    cat(" ", ali_str, "\n", sep="")
}


