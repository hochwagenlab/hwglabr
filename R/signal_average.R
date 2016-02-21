#' Calculate average signal by position
#'
#' Given the ChIP-seq signal on a range of relative genomic positions (e.g. between convergent genes
#' or on ORFs), this function allows you to calculate the average ChIP signal by position. It takes
#' as input a data frame containing the genome-wide signal, for example the output of 'signal_at_conv()'
#' or 'signal_at_orf()'.
#' @param inputData As a data frame containing at least a column named 'position', containing the
#' relative genomic position and a column named 'signal', contianing the corresponding signal.
#' No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current
#' working directory). If 'saveFile = FALSE', output is returned to screen or an R object (if assigned).
#' Defaults to FALSE.
#' @return An R data frame with two columns: position (relative to midpoint of intergenic region) and
#' mean signal.
#' @examples
#' signal_average(WT_conv)
#' 
#' signal_average(WT_S288C_ORFsignal)
#' 
#' signal_average(WT_conv, saveFile = TRUE)
#' @export

signal_average <- function(inputData, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # Make sure the input is a data frame
  if (!is.data.frame(inputData)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'signal_at_conv()' or 'signal_at_orf()' on your data first. Example:\n",
         "WT_signal_dataframe <- signal_at_conv(WT_wiggle)\n",
         "WT_mean_signal <- signal_average(WT_signal_dataframe)", call. = FALSE)
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
                                  "_S288C_mean.txt"), sep = "\t", quote = FALSE)  
    } else {
      write.table(mean_signal, paste0(deparse(substitute(inputData)),
                                  "_SK1_mean.txt"), sep = "\t", quote = FALSE)
    }
  } else {
    return(mean_signal)
  }
  cat(paste0('Completed in ', round((proc.time()[3] - ptm[3]), 1), ' sec.\n'))
}