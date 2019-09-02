#' @export
filter_secondary_alignments <- function(ali, remove_inversions=FALSE){
    if( "tp" %in% names(ali) ){
        if(remove_inversions){
            return(ali[ ali$tp ==  "P"])
        } 
        return(ali[ ali$tp %in% c("I", "i", "P"),])
    }   
    stop("No 'tp' flag in alignment, can filter secondary alignments")
}
