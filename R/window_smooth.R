#' Sliding window smooth
#'
#' This function allows you to smooth data using a sliding window.
#' @param wiggleData As a list of the 16 chr wiggle data (output of readall_tab). No default.
#' @param chrNumber A number representing the chromosome to smooth. No default.
#' @param windowSize A number representing the smoothing window size (number of bp). Defaults to 200.
#' @return An R data frame with two columns: genome position and smoothed signal.
#' @examples
#' window_smooth(WT, 1, 200)
#' 
#' window_smooth(WT, 16, 100)
#' @export

window_smooth <- function(wiggleData, chrNumber, windowSize = 200) {
  # Check reference genome
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  
  check_S288C <- any(grep('chrI.', names(wiggleData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(wiggleData), fixed = TRUE))
  
  if (check_S288C) {
    chrNumber <- paste0('chr', chrom_S288C[chrNumber])
    cat("Smoothing ", chrNumber, ' (mapped to S288C)')
  } else if (check_SK1) {
    chrNumber <- paste0('chr', chrom_SK1[chrNumber])
    cat("Smoothing ", chrNumber, ' (mapped to SK1)')
  } else stop('Did not recognize reference genome.')
  
  listIndex <- grep(paste0(chrNumber, '.'), names(wiggleData), fixed = TRUE)
  chromData <- wiggleData[[listIndex]]
  
  data <- as.data.frame(matrix(0, ncol = 2,
                               nrow = (floor(nrow(chromData)/windowSize)-1)))
  
  for (i in 1:(floor(nrow(chromData)/windowSize) - 1)) {
    data[i, ] <- colMeans(chromData[(windowSize*(i - 1) + 1):(windowSize * i), ])
  }
  
  colnames(data) <- c('position', 'signal')
  return (data)
}