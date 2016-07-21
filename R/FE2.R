#' Microarray probe log2 ratios
#'
#' This function allows you to extract Log2 ratios for all probes in a microarray.
#' It takes as input the raw txt files found in the lab's microarray database
#' (at 'Labshare/HTGenomics/Microarray_database/arrays') or downloaded from GEO.
#' @param fileName A string representing the input file name. No default.
#' @param cy0 A number representing the dye color used for the control sample
#' (\code{3} for cy3 or \code{5} for cy5). No default.
#' @param strain A string representing the output file name (the sample name,
#' typically a yeast strain). No default.
#' @return A flat file with four columns:
#' \enumerate{
#'   \item chr number
#'   \item Probe start position (bp number)
#'   \item Probe end position (bp number)
#'   \item Log2 ratio
#' }
#' Load this file using base R function \code{read.table('/path/to/file', header = TRUE)}.
#' @examples
#' FE2('GSM873122.txt', 3, 'wt_Rec8_1_Cy3')
#' 
#' FE2('247a.txt', 3, '247a')
#' @export

FE2 <- function(fileName, cy0, strain) {
  # Load package "marray"
  if (!requireNamespace("marray", quietly = TRUE)) {
    stop("R package 'marray' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  library(marray)
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using 'data-raw/data_internal.R'; stored in 'R/sysdata.rda'
  #----------------------------------------------------------------------------#
  plat <- SK1rosetta
  maData <- read.Agilent(fnames = fileName, name.Rf = "rMeanSignal",
                         name.Gf = "gMeanSignal", name.Rb = "rBGMeanSignal",
                         name.Gb = "gBGMeanSignal", sep = "\t")
  maNorm <- maNorm(maData, norm = "printTipLoess")
  if (cy0 == 5)
  {
    lr <- -maNorm@maM
  }
  if (cy0 == 3)
  {
    lr <- maNorm@maM
  }
  data <- as.data.frame(matrix(0, nrow = nrow(plat),ncol = 4))
  data[, 1] <- plat[, 'chr']
  data[, 2] <- plat[, 'start']
  data[, 3] <- plat[, 'stop']
  data[, 4] <- lr[, 1]
  data = data[order(data[, 2]), ]
  data = data[order(data[, 1]), ]
  data = data[which(data[, 1] != 0), ]
  data = data[which(data[, 2] != 0), ]
  data = data[which(data[, 3] != 0), ]
  colnames(data) <- c('chr', 'start', 'end', "Log2Ratio")
  write.table(data, strain, col.names = T, row.names = F, quote = F, sep = "\t")
  # save the signal data into a txt file in your working directory.
}