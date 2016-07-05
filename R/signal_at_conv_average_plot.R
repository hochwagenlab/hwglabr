#' Line plot of average signal between convergent genes
#'
#' \strong{Deprecated! Use \code{\link{signal_at_intergen_average_plot}} instead.} \cr
#' This function allows you to make a line plot of the ChIP signal over the average
#' intergenic region centered on midpoints of convergent genes. It takes as input a
#' data frame containing the average signal centered on the midpoints of convergent
#' gene regions.
#' 
#' To generate the input for this function starting from an R list of wiggle data
#' for the 16 chromosomes you should run:
#' 1.  \code{signal_at_conv()} to pull out the signal at every convergent gene region.
#' 
#' 2. \code{signal_average()} to calculate the average signal over all regions.
#' @param inputDataA A data frame of average signal between convergent genes: relative
#' position and average signal. No default.
#' @param inputDataB Optional data in the same format for a second sample. No default.
#' @param genome A string representing the genome used for mapping. This is used in the title
#' of the plot only. No default.
#' @param yMax Optional number to be used as the max Y scale value in the plot.
#' @param onScreen Boolean indicating plots should be returned to the screen
#' (\code{onScreen = TRUE}) or written to .png files (\code{onScreen = FALSE}).
#' Defaults to \code{TRUE}.
#' @param legendXcoord A number representing the X coordinate to locate legend.
#' Defaults to minimum X (left-aligned).
#' @param legendYcoord A number representing the Y coordinate to locate legend.
#' Defaults to maximum Y (top-aligned).
#' @param colorA Optional R color for sample A. Defaults to \code{grey50}.
#' @param colorB Optional R color for sample B. Defaults to \code{green}.
#' @return A line plot of one or two samples, either on screen or as a .png file
#' (in the working directory).
#' @examples
#' signal_at_conv_average_plot(WT_conv_mean_signal, genome = 'S288C')
#' 
#' signal_at_conv_average_plot(WT_conv_mean_signal, dot1_conv_mean_signal, genome = 'SK1',
#'                             yMax = 3, onScreen = FALSE, legendXcoord = -500,
#'                             legendYcoord = 1, colorA = 'red', colorB = 'green')
#' @export

signal_at_conv_average_plot <- function(inputDataA, inputDataB, genome,
                                        yMax, onScreen = TRUE,
                                        legendXcoord = xMin + xMin * 0.2,
                                        legendYcoord = yMax + yMax * 0.05,
                                        colorA = 'grey50', colorB = 'orange') {
  
  cat('Note: This function is deprecated!\n',
      'From version 0.2 on it is no longer maintained. 
      Use signal_at_intergen_average_plot() instead.')
  
  # Make sure the input is the 16x2 data frame returned by chr_cov
  if (!is.data.frame(inputDataA)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'signal_at_conv()' and 'signal_average()' on your data first.
         Example:\n",
         "WT_signal_dataframe <- signal_at_conv(WT_wiggle)\n",
         "WT_mean_signal <- signal_average(WT_signal_dataframe)", call. = FALSE)
  }
  
  if (!missing(inputDataB)) {
    if (!is.data.frame(inputDataB)) {
      stop(deparse(substitute(inputDataB)), " is of wrong format - not an R data frame.\n",
           "Please run 'signal_at_conv()' and 'signal_average()' on your data first.
           Example:\n",
           "WT_signal_dataframe <- signal_at_conv(WT_wiggle)\n",
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
  
  if (missing(yMax)) {
    if (!missing(inputDataB)) {
      # Get highest of the maxima of the two strains to set y axis maximum
      yMax_A <- ceiling(max(dataA[, 2]))
      yMax_B <- ceiling(max(dataB[, 2]))
      yMax <- tail(sort(c(yMax_A, yMax_B)), 1)
      
      # Get lowest of the minima of the two strains to set y axis minimum
      yMin_A <- floor(min(dataA[, 2]))
      yMin_B <- floor(min(dataB[, 2]))
      yMin <- head(sort(c(yMin_A, yMin_B)), 1)
    } else {
      yMax <- ceiling(max(dataA[, 2]))
      yMin <- floor(min(dataA[, 2]))
    }
  }
  
  # Have to use double "[[]]" to extract single values from "tbl_df" dplyr objects
  xMin <- dataA[[1, 1]]
  xMax <- dataA[[nrow(dataA), 1]]
  
  if (!onScreen) {
    png(filename = paste0(deparse(substitute(inputDataA)), "_conv_mean", ".png"),
        width = 700, height = 650, units = 'px')
  }
  
  par(mfrow = c(1, 1), mar = c(11, 12, 8, 4), mgp = c(7, 2, 0))
  plot(0, type = 'l', lwd = 5, xaxt = 'n', yaxt = 'n',
       xlim = c(xMin, xMax), ylim = c(yMin, yMax),
       xlab = "Distance from midpoints of\nintergenic regions (bp)",
       ylab = 'Signal', main = paste0('Convergent genes\n(mapped to ', genome, ' genome)'),
       cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2.5, bty = "n")
  
  axis(1, at = c(xMin, 0, xMax), lab = c(-500, 0, 500),
       las = 1, lwd = 4, cex.axis = 2.5, cex = 3.0)
  axis(2, at = c(yMin, 0, yMax), lwd = 4, las = 2, cex.axis = 2.5, cex = 3.0)
  abline(v = 0, lty = 2, lwd = 2)
  
  lines(dataA, col = colorA, lwd = 8)
  
  if (!missing(inputDataB)) {
    lines(dataB, col = colorB, lwd = 8)
  }
  
  if (missing(inputDataB)) {
    legend(legendXcoord, legendYcoord, deparse(substitute(inputDataA)),
           bty = 'n', cex = 2.5, text.col = colorA)
  } else {
    legend(legendXcoord, legendYcoord,
           c(deparse(substitute(inputDataA)), deparse(substitute(inputDataB))),
           bty = 'n', cex = 2.5, text.col = c(colorA, colorB))
  }
  
  if (!onScreen) {
    dev.off()
  }
}