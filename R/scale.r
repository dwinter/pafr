#' Number formatters for scales in base pairs
#'
#' For use with \code{\link{ggplot2}} 
#' @param x  The data (in base pairs) to be formatted as Gb, Mb or Kb
#' @return  A character vector with scale labels
#' @examples
#' \dontrun{
#' ali <- read_paf(system.file("extdata", "fungi.paf", package="pafr"))
#' doplot(ali) + scale_x_continuous("Genomic position", label=Mb_lab)
#'}
#' @export
#' @rdname Gb_lab
Gb_lab <- function(x) paste(x / 1e9, "Gb")

#' @export
#' @rdname Gb_lab
Mb_lab <- function(x) paste(x / 1e6, "Mb")

#' @export
#' @rdname Gb_lab
Kb_lab <- function(x) paste(x / 1e3, "Kb")
