% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal_at_rDNA.R
\name{signal_at_rDNA}
\alias{signal_at_rDNA}
\title{Signal flanking rDNA}
\usage{
signal_at_rDNA(inputData, saveFile = FALSE)
}
\arguments{
\item{inputData}{As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.}

\item{saveFile}{Boolean indicating whether output should be written to a .txt file (in current working
directory). If \code{saveFile = FALSE}, output is returned to screen or an R object (if assigned).
Defaults to \code{FALSE}.}
}
\value{
A local data frame with four columns:
\enumerate{
  \item \code{chr} Chromosome number
  \item \code{position} Nucleotide coordinate (in normalized total length of 1 kb)
  \item \code{signal} ChIP-seq signal at each position (1 to 1000)
  \item \code{gene} Systematic gene name
}
}
\description{
This function allows you to pull out the ChIP signal flanking the rDNA region on
chromosome 12. The function takes as input the wiggle data as a list of 16 chromosomes.
(output of \code{\link{readall_tab}}).\cr
}
\examples{
\dontrun{
signal_at_rDNA(WT)

signal_at_rDNA(WT, saveFile = TRUE)
}
}
