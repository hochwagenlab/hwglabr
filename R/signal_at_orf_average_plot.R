#' Line plot of average signal on ORFs
#'
#' This function allows you to make a line plot of the ChIP signal over the average
#' ORF. It takes as input a data frame containing the average signal over ORFs
#' (+/- flanking regions) genome-wide.
#' 
#' To generate the input for this function starting from an R list of wiggle data
#' for the 16 chromosomes you should run:
#' 
#' 1. "signal_at_orf()" to pull out the signal at every ORF and flanking regions
#' (1/2 the length of the ORF).
#' 
#' 2. "signal_average()" to calculate the average signal at each relative position
#' over all regions genome-wide.
#' @param inputDataA A data frame of average signal on ORFs +/- 1/2 the length of the ORF:
#' relative position and average signal. No default.
#' @param inputDataB Optional data in the same format for a second sample. No default.
#' @param genome A string representing the genome used for mapping. This is used in the title
#' of the plot only. No default.
#' @param yMin Optional number to be used as the minimum Y scale value in the plot.
#' @param yMax Optional number to be used as the maximum Y scale value in the plot.
#' @param onScreen Boolean indicating plots should be returned to the screen (onScreen = TRUE)
#' or written to .png files (onScreen = FALSE). Defaults to TRUE.
#' @param legend_Xcoord A number representing the X coordinate to locate legend.
#' Defaults to minimum X (left-aligned).
#' @param legend_Ycoord A number representing the Y coordinate to locate legend.
#' Defaults to maximum Y (top-aligned).
#' @param colorA Optional R color for sample A. Defaults to 'grey50'.
#' @param colorB Optional R color for sample B. Defaults to 'green'.
#' @return A dot plot of one or two samples, either on screen or as a png file
#' (in the working directory).
#' @examples
#' signal_at_orf_average_plot(WT_orf_S288C_mean_signal, genome = 'S288C')
#' 
#' signal_at_orf_average_plot(WT_orf_mean_signal, dot1_orf_mean_signal, genome = 'SK1',
#'                            yMax = 3, onScreen = FALSE, legend_Xcoord = -500,
#'                            legend_Ycoord = 1, colorA = 'red', colorB = 'green')
#' @export

signal_at_orf_average_plot <- function(inputDataA, inputDataB, genome,
                                       yMin, yMax, onScreen = TRUE,
                                       legend_Xcoord = -200,
                                       legend_Ycoord = yMax + yMax * 0.05,
                                       colorA = 'grey50', colorB = 'orange') {
  
  # Make sure the input is the 16x2 data frame returned by chr_cov
  if (!is.data.frame(inputDataA)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'signal_at_orf()' and 'signal_average()' on your data first.
         Example:\n",
         "WT_signal_dataframe <- signal_at_orf(WT_wiggle)\n",
         "WT_mean_signal <- signal_average(WT_signal_dataframe)", call. = FALSE)
  }
  
  if (!missing(inputDataB)) {
    if (!is.data.frame(inputDataB)) {
      stop(deparse(substitute(inputDataB)), " is of wrong format - not an R data frame.\n",
           "Please run 'signal_at_orf()' and 'signal_average()' on your data first.
           Example:\n",
           "WT_signal_dataframe <- signal_at_orf(WT_wiggle)\n",
           "WT_mean_signal <- signal_average(WT_signal_dataframe)", call. = FALSE)
    }
  }
    
  # Not changing variable names here causes problems with deparse(substitute())
  # in plot legend
  dataA <- inputDataA
  if (!missing(inputDataB)) {
    dataB <- inputDataB
  }
  
  ### Plot
  cat('Plotting... \n')

  if (missing(yMin)) {
    if (!missing(inputDataB)) {
      # Get highest of the maxima of the two strains to set y axis maximum
      yMin_A <- min(dataA[, 2])
      yMin_B <- min(dataB[, 2])
      yMin <- tail(sort(c(yMin_A, yMin_B)), 1)
    } else yMin <- min(dataA[, 2])
  }
    
  if (missing(yMax)) {
    if (!missing(inputDataB)) {
      # Get highest of the maxima of the two strains to set y axis maximum
      yMax_A <- ceiling(max(dataA[, 2]))
      yMax_B <- ceiling(max(dataB[, 2]))
      yMax <- tail(sort(c(yMax_A, yMax_B)), 1)
    } else yMax <- ceiling(max(dataA[, 2]))
  }
  
  if (!onScreen) {
    png(filename = paste0(deparse(substitute(inputDataA)), "_conv_mean", ".png"),
        width = 700, height = 650, unit = 'px')
  }
  
  par(mfrow = c(1, 1), mar = c(11, 12, 8, 4), mgp = c(7, 2, 0))
  plot(0, type = 'l', lwd = 5, xaxt = 'n', yaxt = 'n',
       xlim = c(0, 2000), ylim = c(yMin, yMax),
       xlab = "Normalized position",
       ylab = 'Signal', main = paste0('ORFs (mapped to ', genome, ' genome)'),
       cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2.5, bty = "n")
  
  axis(1, at = c(0, 2000), lab = c('', ''),
       las = 1, lwd = 4, cex.axis = 2.5, cex = 3.0)
  axis(2, at = c(0, yMax), lwd = 4, las = 2, cex.axis = 2.5, cex = 3.0)
  abline(v = c(500, 1500) , lty= 2, lwd = 2)
  
  lines(dataA, col = colorA, lwd = 8)
  
  if (!missing(inputDataB)) {
    lines(dataB, col = colorB, lwd = 8)
  }
  
  if (missing(inputDataB)) {
    legend(legend_Xcoord, legend_Ycoord, deparse(substitute(inputDataA)),
           bty = 'n', cex = 2.5, text.col = colorA)
  } else {
    legend(legend_Xcoord, legend_Ycoord,
           c(deparse(substitute(inputDataA)), deparse(substitute(inputDataB))),
           bty = 'n', cex = 2.5, text.col = c(colorA, colorB))
  }
  
  if (!onScreen) {
    dev.off()
  }
}