#' @export
filter_secondary_alignments <- function(ali){
    if( "tp" %in% names(ali) ){
        return(ali[ ali$tp == "P",])
    }   
    stop("No 'tp' flag in alignment, can filter secondary alignments")
}
