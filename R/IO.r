#' Read genomic intervals in bed format
#'
#' The first three columns of the file specified by \code{file_name} must
#' contain data in the standard bed format (i.e., a genomic interval
#' represented by 0-based half-open interval with seq-id, start and end position). 
#' These columns will be renamed to 'chrom', 'start' and 'end', respectively. Any
#' other columns present in the data will be left unmodified.
#' 
#' The file is read into memory with \code{\link{read.table}}, with the
#' argument \code{sep} set to \code{'\t'} and \code{stringsAsFactors} set to
#' FALSE. All other arguments are left as default, but arguments can be passed
#' from \code{read_bed} to \code{read.table}.

#' @param file_name  Path to the bed file to be read in
#' @param tibble logical  If TRUE, the genomic intervals are returned as
#' a tidy \code{tbl_df}.
#' @param ...  Other arguments passed to \code{read.table}
#' @return Either a \code{data.frame} or a \code{tbl_df} with at least three
#' columns named 'chrom', 'start' and 'end'
#' @importFrom tibble as_tibble
#' @importFrom utils read.table
#' @examples
#' bed_path <- system.file("extdata", "Q_centro.bed", package="pafr")
#' centro <- read_bed(bed_path)
#' centro
#' # Can pass arguments to read.table
#' miss_two <- read_bed(bed_path, skip=2)
#' miss_two 
#' @export
read_bed <- function(file_name, tibble = FALSE, ...) {
    res <- read.table(file_name, sep = "\t", stringsAsFactors = FALSE, ...)
    if (ncol(res) < 3) {
        stop("Bed file must have at least three columns, data has", ncol(res))
    }
    names(res)[1:3] <- c("chrom", "start", "end")
    if (tibble) {
        return(as_tibble(res))
    }
    res
}


.make_numeric <- function(df, cols) {
    for (i in cols) {
        df[[i]] <- as.numeric(df[[i]])
    }
    df
}



#' @importFrom dplyr bind_rows
process_tags <- function(raw_tags){
    tags <- character()
    to_numeric <- integer()
    res <- list()
    n <- length(raw_tags)
    t_idx <- 0
    for(ali_idx in seq_along(raw_tags)) {
        split_tags <- str_split(raw_tags[[ali_idx]], ":")
        for(tag in split_tags){
            if ( !(tag[1] %in% tags) ){
                t_idx <- t_idx + 1
                if( tag[2] %in% c("f", "H", "i")){
                    to_numeric <- c(to_numeric, t_idx)
                }
                res[[ tag[1] ]] <- rep(NA, n)
                tags <- c(tags, tag[1])
            }
            res[[ tag[1] ]][ali_idx] <- tag[3]
        }
    }
    
    for(i in to_numeric){
        res[[i]] <- as.numeric(res[[i]])
    }
    bind_rows(res)
}

#' Read a genomic alignment in PAF format
#'
#' See the package vignette for detailed information on the file format and its
#' representation as an R object.
#
#' @param file_name  Path to the .paf file
#' @param tibble logical  If TRUE, the genomic alignments are returned as
#' a tidy \code{tbl_df}
#' @param include_tags  logical if TRUE (default) read additional information
#' about each alignment encoded as PAF tags. Setting this to FALSE will speed up
#' parsing of paf alignments, specially those with large CIGAR strings/
#' @return Either a \code{pafr} object, which acts as a \code{data.frame}, or a
#' \code{tbl_df} containing information on genomic alignments. The contents of
#' this table are described in detail in the pafr package vingette.
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @examples
#' ali <- read_paf( system.file("extdata", "fungi.paf", package="pafr") )
#' ali
#' @export
read_paf <- function(file_name, tibble=FALSE, include_tags=TRUE) {
    lines <- readLines(file_name)
    if( include_tags ){
        tokens <- str_split(lines, "\t")

    } else {
        tokens <-  str_split(lines, "\t", 13)
    }
    paf_cols <- lapply(tokens, "[", 1:12)
    col_names <- c("qname", "qlen", "qstart", "qend", "strand",
                   "tname", "tlen", "tstart", "tend", "nmatch",
                   "alen", "mapq")
    #bind_rows is a decent amount faster then do.call(rbind.data.frame, ...),
    # but requres a named vector.. because reasons? Even with this extra set
    # this is faster than rbind.data.frame
    for(i in seq_along(paf_cols)){
        attr(paf_cols[[i]], "names") <- col_names
    }
    res <- bind_rows(paf_cols)
    res <- .make_numeric(as.data.frame(res), c(2, 3, 4, 7:12))
    #first 12 columns are always the paf alignment. Anything after this is a
    #tag.
    if(include_tags){
      if (any(lengths(tokens) > 12)) {
            raw_tags <- sapply(tokens, function(x) paste(x[13:length(x)]))
            res <- cbind.data.frame(res, process_tags(raw_tags))
        }
    }
    if (tibble) {
        return(as_tibble(res))
    }
    class(res) <- c("pafr", class(res))
    res
}

#' Coerce a data.frame or tibble into a pafr object
#' 
#' The main reason to use this function is speed up the process of reading in
#' a large paf file that has no tags. Functions like read.table, read_delim
#' (reader) and fread (data.table) can process a 12 column file  more quickly
#' than pafr's read_paf. If you you do not need tag data for your analyses or
#' visualizations, it might make sense to use a fast reading function to get a
#' 12 column data.frame,  convert that data.frame into a `pafr  object with this
#' function. The `pafr` object can then work easily with the functions in this package.
#' @export 
#' @param paf_data_frame a data.frame object with 12 columns. Column names and
#' types will be over-written in the returned object
#' @return a pafr object 
#' @seealso read_paf
as_paf <- function(paf_data_frame){
    #lets be assertive
    if( ncol(paf_data_frame) != 12 ){
        stop("data.frame should be in paf format, with 12 columns")
    }
    col_names <- c("qname", "qlen", "qstart", "qend", "strand",
                   "tname", "tlen", "tstart", "tend", "nmatch",
                   "alen", "mapq")
    for(col_idx in c(2:4, 7:12)){
        if( !(is_numericable(paf_data_frame[[col_idx]]) ) ){
            msg <- paste0("Column", col_idx, "('", col_names[col_idx], 
                          "should be numeric or convertable to numeric without creating NAs")
            stop(msg)
        }
    }
    #ok, looks good
    res <- .make_numeric(paf_data_frame, c(2, 3, 4, 7:12))
    names(res) <- col_names
    class(res) <- c("pafr", "data.frame")
    res
}    

#is a vector, possibly in character form, able to be treated as a numeric?
is_numericable <- function(vec){
    UseMethod("is_numericable", vec)
}

is_numericable.character <- function(vec){
    suppressWarnings(res <- as.numeric(vec))
    !(any(is.na(res)))
}

is_numericable.numeric <- function(vec){
    TRUE
}

is_numericable.factor <- function(vec){
    is_numericable(as.character(vec))
}
#' @export
print.pafr <- function(x, ...) {
    if (nrow(x) == 1) {
        print_alignment(x)
        return(invisible(x))
    }
    print_pafr_rows(x)
    return(invisible(x))
}
    
print_pafr_rows <- function(x) {
    cat("pafr object with ", nrow(x), " alignments (",
        round(sum(x$alen) / 1e6, 1), "Mb)\n", sep = "")
    cat(" ", length(unique(x$qname)), " query seqs\n", sep = "")
    cat(" ", length(unique(x$tname)), " target seqs\n", sep = "")
    ntags <- ncol(x) - 12
    if (ntags < 12) {
        if (ntags > 0) {
            tags <- paste(names(x)[13:ncol(x)], collapse = ", ")
        } else {
            tags <- ""
        }

    } else {
        tags <- paste(paste(names(x)[13:23], collapse = ", "), "...")
    }
    cat(" ", ntags, " tags: ", tags, "\n", sep = "")
}

print_alignment <- function(x) {
    ali_str <- paste0(x[["qname"]], ":", x[["qstart"]], "-", x[["qend"]], " v ",
                      x[["tname"]], ":", x[["tstart"]], "-", x[["tend"]])
    cat("Single pafr alignment:\n")
    cat(" ", ali_str, "\n", sep = "")
}
