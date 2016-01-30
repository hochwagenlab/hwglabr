#' Calculate average signal between convergent genes
#'
#' This function allows you to calculate the ChIP signal over the average intergenic region centered
#' on midpoints of convergent genes. It takes as input a data frame containing the genomewide signal
#' centered on the midpoints of convergent gene regions (output of 'signal_at_conv()') and calculates
#' the average for each relative position.
#' @param inputData As a data frame (output of 'signal_at_conv()'). No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current
#' working directory). If 'saveFile = FALSE', output is returned to screen or an R object (if assigned).
#' Defaults to FALSE.
#' @return An R data frame with two columns: genome position and mean signal.
#' @examples
#' average_signal_at_conv(WT_conv)
#' 
#' average_signal_at_conv(WT_conv, region_size = 1500, saveFile = TRUE, input_dataFrame = FALSE)
#' @export

average_signal_at_conv <- function(inputData, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # Make sure the input is a data frame
  if (!is.data.frame(inputData)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'signal_at_conv()' on your data first. Example:\n",
         "WT_signal_dataframe <- signal_at_conv(WT_wiggle)\n",
         "average_signal_at_conv(WT_signal_dataframe)", call. = FALSE)
    }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  # Calculate averages for each relative position
  mean_signal <- inputData %>%
    group_by(position) %>%
    summarise(Mean_signal = mean(signal))
  
  if(saveFile) {
    cat(paste0('Saving file...\n'))
    if(check_S288C) {
      write.table(mean_signal, paste0(deparse(substitute(inputData)),
                                  "_S288C_conv_mean.txt"), sep = "\t", quote = FALSE)  
    } else {
      write.table(mean_signal, paste0(deparse(substitute(inputData)),
                                  "_SK1_conv_mean.txt"), sep = "\t", quote = FALSE)
    }
  } else {
    return(mean_signal)
  }
  cat(paste0('Completed in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.\n'))
}