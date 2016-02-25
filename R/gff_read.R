#' Read gff file
#'
#' This function allows you to read a gff file. It uses \code{read.table()} to load the
#' gff as a data frame, adds the column names and sets the correct column data classes.
#' Adapted from function in the \href{https://bioconductor.org/packages/release/data/experiment/html/davidTiling.html}{davidTiling experimental package}.
#' Also posted on the Bioconductor support system \href{https://support.bioconductor.org/p/24657/}{here}.
#' @param gffFile The path to the Gff file to read. No default.
#' @param nRows Optional maximum number of rows to read (the exact argument to \code{nRows}
#' of \code{read.table()}). Defaults to all rows.
#' @return Gff as data frame with appropriate column names and column data set to the
#' correct class.
#' @examples
#' gff_read(s288C_annotation_R64_modified.gff)
#' @export

gff_read <- function(gffFile, nRows = -1) {
  cat("Reading ", gffFile, ": ", sep = "")
  gff = read.table(gffFile, sep = "\t", as.is = TRUE, quote = "",
                   header = FALSE, comment.char = "#", nrows = nRows,
                   colClasses = c("character", "character", "character",
                                  "integer", "integer", "character", "character",
                                  "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}
