% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signal_at_orf_se.R
\name{signal_at_orf_se}
\alias{signal_at_orf_se}
\title{Signal around ORF start or end positions}
\usage{
signal_at_orf_se(inputData, gff, gffFile, limit = "start", upstrExt = 500,
  downstrExt = 1500, saveFile = FALSE)
}
\arguments{
\item{inputData}{As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.}

\item{gff}{Optional dataframe of the gff providing the ORF cordinates. Must be provided if
\code{gffFile} is not. No default. Note: You can use the function \code{\link{gff_read}} in hwglabr to
load your selected gff file.}

\item{gffFile}{Optional string indicating path to the gff file providing the ORF cordinates. Must be
provided if \code{gff} is not. No default.}

\item{limit}{String indicating whether to use translation start or end as a the
reference point. Accepts one of \code{start} or \code{end}. Defaults to \code{start}.}

\item{upstrExt}{Number specifying the extension in bp to collect upstream of the
reference limit. Defaults to \code{500}.}

\item{downstrExt}{Number specifying the extension in bp to collect downstream of the
reference limit. Defaults to \code{1500}.}

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
This function allows you to pull out the ChIP signal around the translation start
or end for all ORFs in the genome. It collects the signal between specified extensions
up and downstream of the chosen ORF limit.\cr
The function takes as input the wiggle data as a list of 16 chromosomes.
(output of \code{\link{readall_tab}}).
\cr \cr
\strong{Note:} Our wiggle data always contains gaps with missing chromosome coordinates
and ChIP-seq signal. The way this function deals with that is by skipping affected genes.
The number of skipped genes in each chromosome is printed to the console, as well as the
final count (and percentage) of skipped genes. \cr
}
\examples{
\dontrun{
signal_at_orf_se(WT, gff = gff)

signal_at_orf_se(WT, gffFile = S288C_annotation_modified.gff, limit = 'end',
                 upstrExt = 1500, downstrExt = 500, saveFile = TRUE)
}
}
