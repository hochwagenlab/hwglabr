#' Compute mean signal and SD by position of output of \code{\link{signal_from_telomeres}}
#'
#' Given the  output of \code{\link{signal_from_telomeres}}, this function allows
#' you to compute the mean ChIP-seq signal and  standard deviation by position (as
#' distance to telomere) either for all chromosomes together or in two separate groups
#' for small (chrs 1, 3, and 6) and large chromosomes (all remaining chromosomes).\cr\cr
#' @param inputData A list produced by calling \code{\link{signal_from_telomeres}}. No default.
#' @param separateSmallAndLarge Boolean specifying whether to generate a single data
#' frame for all chromosomes (\code{separateSmallAndLarge = FALSE}) or two data frames
#' separating small (chrs 1, 3, and 6) and large chromosomes (all remaining chromosomes).
#' Defaults to \code{FALSE}.
#' @return If \code{separateSmallAndLarge = FALSE}, a dataframe with three columns:
#' \enumerate{
#'   \item \code{distance_to_telomere} Distance to telomere in bp
#'   \item \code{mean_signal} Average ChIP-seq signal
#'   \item \code{sd} Standard deviation
#' }
#' If \code{separateSmallAndLarge = TRUE}, a list of two dataframes (each with the
#' same three columns explained above):
#' \enumerate{
#'   \item \code{small_chrs} Data frame of chromosomes 1, 3 and 6
#'   \item \code{large_chrs} Data frame of remaining chromosomes
#' }
#' @examples
#' \dontrun{
#' signal_from_telomeres_average(WT_signal)
#' 
#' signal_from_telomeres_average(signalFromTelomeres, separateSmallAndLarge = TRUE)
#' }
#' @export

signal_from_telomeres_average <- function(inputData, separateSmallAndLarge = FALSE){
  ptm  <- proc.time()
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  message('"rbind"ing data frames...')
  # Make data frame from lists of chromosome arms
  inputData$small_chrs <- do.call('rbind', inputData$small_chrs)
  inputData$large_chrs <- do.call('rbind', inputData$large_chrs)
  
  if(separateSmallAndLarge){
    message('\nCalculating mean and SD...')
    # Group by genome position and average signal
    inputData$small_chrs <- dplyr::summarise(dplyr::group_by(inputData$small_chrs,
                                                             distance_to_telomere),
                                             mean_signal = mean(signal, na.rm = TRUE),
                                             sd = sd(signal, na.rm = TRUE))
    inputData$large_chrs <- dplyr::summarise(dplyr::group_by(inputData$large_chrs,
                                                             distance_to_telomere),
                                             mean_signal = mean(signal, na.rm = TRUE),
                                             sd = sd(signal, na.rm = TRUE))
    finalData <- inputData
    
  } else {
    finalData <- do.call('rbind', inputData)
    
    # Group by genome position and average signal
    message('\nCalculating mean and SD...')
    finalData <- dplyr::summarise(dplyr::group_by(finalData, distance_to_telomere),
                                  mean_signal = mean(signal, na.rm = TRUE),
                                  sd = sd(signal, na.rm = TRUE))
  }
  
  message(paste0('\n\nCompleted in ', round((proc.time()[3] - ptm[3]), 2), ' sec.'))
  return(finalData)
}