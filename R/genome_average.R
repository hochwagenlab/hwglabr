#' Calculate genome-wide average signal
#'
#' Given the wiggle data as a list of 16 chromosomes, this function returns a single
#' value corresponding to the genome-wide signal average.
#' @param inputData List of the 16-chr wiggle data (output of \code{\link{readall_tab}}).
#' No default.
#' @return A numeric value corresponding to the genome-wide average signal.
#' @examples
#' \dontrun{
#' average_signal <- genome_average(WT)
#' }
#' @export

genome_average <- function(inputData) {
  # Make sure the input is a list with 16 elements
  if (!is.list(inputData) | length(inputData) != 16) {
    stop("Wrong input data - not a list with 16 elements.\n",
         "Please provide a list of 16 data frames, one for each of the 16 chromosomes",
         "(the output of readall_tab).", call. = FALSE)
    }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  library(dplyr)
  message('Genome-wide signal average for sample ',
          deparse(substitute(inputData)), ': ',
          appendLF = FALSE)
  #average <- mean(unlist(sapply(inputData, function(x) x[, 2])), na.rm=TRUE)
  average <- colMeans(do.call("bind_rows", inputData))[2]
  message(round(average, 3))
  
  return(as.numeric(average))
}