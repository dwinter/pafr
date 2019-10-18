% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_ali.r
\name{filter_secondary_alignments}
\alias{filter_secondary_alignments}
\title{Remove secondary alignments form a pafr alignment}
\usage{
filter_secondary_alignments(ali, remove_inversions = FALSE)
}
\arguments{
\item{ali}{genomic alignment in \code{pafr} or \code{tbl_df} format, as
returned by \code{\link{read_paf}}}

\item{remove_inversions, }{logical, if \code{TRUE} also remove inversions (tp
flag 'I' or 'i') from alignment.}
}
\description{
Remove secondary alignments form a pafr alignment
}
\examples{
ali <- read_paf( system.file("extdata", "fungi.paf", package = "pafr") ) 
ali
filter_secondary_alignments(ali)
}