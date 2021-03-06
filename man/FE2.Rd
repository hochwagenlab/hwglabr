% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FE2.R
\name{FE2}
\alias{FE2}
\title{Microarray probe log2 ratios}
\usage{
FE2(fileName, cy0, strain)
}
\arguments{
\item{fileName}{A string representing the input file name. No default.}

\item{cy0}{A number representing the dye color used for the control sample
(\code{3} for cy3 or \code{5} for cy5). No default.}

\item{strain}{A string representing the output file name (the sample name,
typically a yeast strain). No default.}
}
\value{
A flat file with four columns:
\enumerate{
  \item chr number
  \item Probe start position (bp number)
  \item Probe end position (bp number)
  \item Log2 ratio
}
Load this file using base R function \code{read.table('/path/to/file', header = TRUE)}.
}
\description{
This function allows you to extract Log2 ratios for all probes in a microarray.
It takes as input the raw txt files found in the lab's microarray database
(at 'Labshare/HTGenomics/Microarray_database/arrays') or downloaded from GEO.
}
\examples{
\dontrun{
FE2('GSM873122.txt', 3, 'wt_Rec8_1_Cy3')

FE2('247a.txt', 3, '247a')
}
}
