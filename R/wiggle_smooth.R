#' Sliding window smooth of wiggle data
#'
#' This function allows you to smooth wiggle data using a sliding window.
#' @param wiggleData As a list of the 16 chr wiggle data (output of \code{readall_tab()}). No default.
#' @param chrNumber A number representing the chromosome to smooth. No default.
#' @param bandwidth A number representing the length of the smoothing window in bp
#' (or the Gaussian kernel bandwith, if \code{useKsmooth = TRUE}). Defaults to 200.
#' @param useKsmooth Boolean indicating choice of smoothing function:
#' \enumerate{
#'   \item \code{useKsmooth = FALSE}: use a simple sliding window smoother. Smoothing is performed by
#'   sliding a window of the specified size (\code{bandwidth} argument) over all genomic positions in
#'   the data and replacing the position values by the middle position and the signal values by
#'   their mean.
#'   \item \code{useKsmooth = TRUE}: use a Gaussian Kernel Regression Smoother. Smoothing is performed
#'   using function \code{ksmooth()} from 'stats' package using the default \code{normal} kernel and the
#'   specified bandwith.
#'   }
#'   Defaults to \code{FALSE}.
#' @return An R data frame with two columns: genome position and smoothed signal.
#' @examples
#' wiggle_smooth(WT, 1, 200)
#' 
#' wiggle_smooth(WT, 16, 100)
#' @export

wiggle_smooth <- function(wiggleData, chrNumber, bandwidth = 200, useKsmooth = FALSE) {
  # Check reference genome
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  
  check_S288C <- any(grep('chrI.', names(wiggleData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(wiggleData), fixed = TRUE))
  
  if (check_S288C) {
    chrNumber <- paste0('chr', chrom_S288C[chrNumber])
    cat("Smoothing ", chrNumber, ' - mapped to S288C\n(Chrs numbered using roman numerals)')
  } else if (check_SK1) {
    chrNumber <- paste0('chr', chrom_SK1[chrNumber])
    cat("Smoothing ", chrNumber, ' - mapped to SK1\n(Chrs numbered using arabic numerals)')
  } else stop('Did not recognize reference genome.')
  
  listIndex <- grep(paste0(chrNumber, '.'), names(wiggleData), fixed = TRUE)
  chromData <- wiggleData[[listIndex]]
  
  if (!useKsmooth) {
    data <- as.data.frame(matrix(0, ncol = 2,
                                 nrow = (floor(nrow(chromData)/bandwidth)-1)))
    
    for (i in 1:(floor(nrow(chromData)/bandwidth) - 1)) {
      data[i, ] <- colMeans(chromData[(bandwidth*(i - 1) + 1):(bandwidth * i), ])
    }
  } else {
    dataList <- ksmooth(as.data.frame(chromData)[, 1], as.data.frame(chromData)[, 2],
                        bandwidth = bandwidth)
    data <- cbind(dataList[[1]], dataList[[2]])
  }
  
  
  colnames(data) <- c('position', 'signal')
  return (data.frame(data))
}
