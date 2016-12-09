#' Calculate genome-wide signal average in wiggle data
#'
#' Given the ChIP-seq wiggle data as a list of 16 chromosomes (output of
#' \code{\link{readall_tab}}), this function calculates the genome-wide signal average.
#' @param inputData As a list of the 16 chr wiggle data (output of
#' \code{\link{readall_tab}}). No default.
#' @return A decimal value of type numeric.
#' @examples
#' \dontrun{
#' wiggle_signal_average(wt_wiggle)
#' }
#' @export

wiggle_signal_average <- function(inputData) {
  # Make sure the input is a list with 16 elements
  if (!is.list(inputData) | !length(inputData) == 16) {
    stop("Wrong input data - not a list with 16 elements.\n",
         "The data should be the output of 'readall_tab()'. Example:\n",
         "WT_wiggle <- readall_tab(readall_tab('/Path/to/WT/wiggle/folder/'))\n",
         "WT_mean_signal <- wiggle_signal_average(WT_wiggle)", call. = FALSE)
  }
  
  # Check genome just to let user know; not used in calculation
  check_S288C <- any(grep("chrI.", names(inputData), fixed = TRUE))
  check_SK1 <- any(grep("chr01.", names(inputData), fixed = TRUE))
  
  # Check reference genome and load respective chromosome number vector
  if (check_S288C) {
    message('Detected ref. genome - S288C')
  } else if (check_SK1) {
    print('Detected ref. genome - SK1')
  } else stop("Did not recognize reference genome.
              Check that chromosome numbers are in the usual format, e.g. 'chrI' or 'chr01'.")
  
  mean_signal <- sum(sapply(inputData, function(x) sum(x[, 2]))) / sum(sapply(inputData, nrow))
  message('Calculating genome-wide average signal')
  message('Done!')
  
  return(mean_signal)
}
