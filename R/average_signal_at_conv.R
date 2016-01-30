#' Calculate average signal between convergent genes
#'
#' This function allows you to calculate the ChIP signal over the average intergenic region centered
#' on midpoints of convergent genes. It takes as input a data frame containing the genome-wide signal
#' centered on the midpoints of convergent gene regions (output of 'signal_at_conv()') and calculates
#' the average for each relative position (e.g. between midpoint of intergenic convergent gene regions
#' + and - 500 bp).
#' @param inputData As a data frame (output of 'signal_at_conv()'). No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current
#' working directory). If 'saveFile = FALSE', output is returned to screen or an R object (if assigned).
#' Defaults to FALSE.
#' @return An R data frame with two columns: position (relative to midpoint of intergenic region) and
#' mean signal.
#' @examples
#' average_signal_at_conv(WT_conv)
#' 
#' average_signal_at_conv(WT_conv, saveFile = TRUE)
#' @export

average_signal_at_conv <- function(inputData, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # Make sure the input is a data frame
  if (!is.data.frame(inputData)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'signal_at_conv()' on your data first. Example:\n",
         "WT_signal_dataframe <- signal_at_conv(WT_wiggle)\n",
         "WT_mean_signal <- average_signal_at_conv(WT_signal_dataframe)", call. = FALSE)
    }
  
  # Check reference genome and load appropriate convergent gene regions 
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  check_S288C <- any(inputData[, 1] == 'chrI')
  check_SK1 <- any(grep('chr01.', names(inputData), fixed = TRUE))
  
  if (check_S288C) {
    cat('Detected ref. genome - S288C\n')
    # Load the data:
    data("S288C_conv")
    conv <- S288C_conv
    chrom <- chrom_S288C
    
  } else if (check_SK1) {
    print('Detected ref. genome - SK1')
    # Load the data:
    data("SK1_conv")
    conv <- SK1_conv
    chrom <- chrom_SK1
  } else stop("Did not recognize reference genome.
               Check that chromosome numbers are in the usual format, e.g. 'chrI' or 'chr01'.")
  
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