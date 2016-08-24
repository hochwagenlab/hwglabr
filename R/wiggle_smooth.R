#' Sliding window smooth of wiggle data
#'
#' This function allows you to smooth wiggle data using a sliding window.
#' @param wiggleData Accepts input in one of the following formats:
#' \enumerate{
#'   \item An R list of the 16 chromosome wiggle data (output of \code{readall_tab()}).
#'   \item An element of an R list of the form described in the first item above.
#'   Can be extracted with either '[]' or '[[]]'.
#'   \item An R data frame in the same format as the individual chromosome data
#'   frames composing the list described in the first item above.
#'   }
#'   No default.
#' @param chrNumber An integer representing the chromosome to smooth. Will be ignored
#' in case the provided \code{wiggleData} is not a list of data for the 16 chromosomes.
#' No default.
#' @param bandwidth An integer representing the length of the smoothing window in bp
#' (or the Gaussian kernel bandwith, if \code{useKsmooth = TRUE}). Defaults to \code{200}.
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
#' wiggle_smooth(rec8, 16, 100)
#' 
#' wiggle_smooth(WT[[5]], bandwidth = 1000)
#' 
#' wiggle_smooth(WT[[9]], bandwidth = 1000, useKsmooth = T)
#' @export

wiggle_smooth <- function(wiggleData, chrNumber, bandwidth = 200, useKsmooth = FALSE) {
  # List of chr numbers for both S288c and SK1
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  
  # Get data, while checking if input is list of chromosomes or one chromosome table
  list_of_16 <- is.list(wiggleData) & !is.data.frame(wiggleData) & length(wiggleData) == 16
  element_from_list_of_16 <- is.list(wiggleData) & !is.data.frame(wiggleData) & length(wiggleData) == 1
  df_from_list_of_16 <- is.list(wiggleData) & is.data.frame(wiggleData)
  if(list_of_16){ # Was a list of 16 chrs provided?
    if(any(grep('chrI.', names(wiggleData), fixed = TRUE))){
      listIndex <- grep(paste0('chr', chrom_S288C[chrNumber], '.'),
                        names(wiggleData), fixed = TRUE)
      chromData <- wiggleData[[listIndex]] 
    } else if(any(grep('chr01.', names(wiggleData), fixed = TRUE))){
      listIndex <- grep(paste0('chr', chrom_SK1[chrNumber], '.'),
                        names(wiggleData), fixed = TRUE)
      chromData <- wiggleData[[listIndex]] 
    }
  } else if(element_from_list_of_16){ # Was a selected element from the list (e.g. [1]) provided?
    chromData <- wiggleData[[1]]
    if(!missing(chrNumber)){
      message('\nNote: You provided a chromosome number, but data for a specific
          chromosome only. The chromosome number (', chrNumber, ') will be ignored!\n')
    }
    } else if(df_from_list_of_16 | !is.list(wiggleData)){
      # Was a selected df extracted from the list (e.g. [[1]]) or a df provided?
      chromData <- wiggleData
      if(!missing(chrNumber)){
        message('\nNote: You provided a chromosome number, but data for a specific
            chromosome only. The chromosome number (', chrNumber, ') will be ignored!\n')
      }
      } else stop('Could not recognize data. Please check that you provided either a list
of 16 data frames (one per chromosome) and a "chrNumber" or a list element or data frame
for a selected chromosome.')
  
  message("\nSmoothing data...")
  
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
  message("\nDone!")
  return (data.frame(data))
}