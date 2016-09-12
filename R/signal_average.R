#' Calculate average signal by position
#'
#' Given the ChIP-seq signal on a range of relative genomic positions (e.g. at intergenic regions
#' or on ORFs), this function allows you to calculate the average ChIP signal by position. It takes
#' as input a data frame containing the genome-wide signal, for example the output of
#' \code{\link{signal_at_intergen}} or \code{\link{signal_at_orf}}.
#' @param inputData As a data frame containing at least a column named \code{position}, containing the
#' relative genomic position and a column named \code{signal}, containing the corresponding signal.
#' No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current
#' working directory). If \code{saveFile = FALSE}, output is returned to screen or an R object
#' (if assigned). Defaults to \code{FALSE}.
#' @return An R data frame with two columns: position (relative genome coordinate) and
#' mean_signal (average signal at each relative coordinate).
#' @examples
#' \dontrun{
#' signal_average(WT_conv)
#' 
#' signal_average(WT_S288C_ORFsignal)
#' 
#' signal_average(WT_conv, saveFile = TRUE)
#' }
#' @export

signal_average <- function(inputData, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # Make sure the input is a data frame
  if (!is.data.frame(inputData)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'signal_at_intergen()' or 'signal_at_orf()' on your data first. Example:\n",
         "WT_signal_dataframe <- signal_at_intergen(WT_wiggle)\n",
         "WT_mean_signal <- signal_average(WT_signal_dataframe)", call. = FALSE)
    }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  check_S288C <- any(inputData[, 1] == 'chrI' | inputData[, 1] == 'chrII')
  check_SK1 <- any(inputData[, 1] == 'chr01' | inputData[, 1] == 'chr02')
  
  # Check reference genome and load respective chromosome number vector
  if (check_S288C) {
    message('Detected ref. genome - S288C\n')
  } else if (check_SK1) {
    print('Detected ref. genome - SK1')
  } else stop("Did not recognize reference genome.
              Check that chromosome numbers are in the usual format, e.g. 'chrI' or 'chr01'.")
  
  # Calculate averages for each relative position
  #mean_signal <- inputData %>%
  #  group_by(position) %>%
  #  summarise(mean_signal = mean(signal))
  
  mean_signal <- dplyr::summarise(dplyr::group_by(inputData, position),
                                  mean_signal = mean(signal, na.rm = TRUE))
  
  if(saveFile) {
    message(paste0('Saving file...'))
    if(check_S288C) {
      write.table(mean_signal, paste0(deparse(substitute(inputData)),
                                  "_S288C_mean.txt"), sep = "\t",
                  quote = FALSE, row.names = FALSE)
    } else {
      write.table(mean_signal, paste0(deparse(substitute(inputData)),
                                  "_SK1_mean.txt"), sep = "\t",
                  quote = FALSE, row.names = FALSE)
    }
  } else {
    return(mean_signal)
  }
  message(paste0('Completed in ', round((proc.time()[3] - ptm[3]), 1), ' sec.'))
}
