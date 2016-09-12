#' Wiggle data line plot
#'
#' This function allows you to plot wiggle data for a selected chromosome.
#' It takes as input either the raw wiggle data (for example one element of the 16 chromosome list)
#' or the output of \code{wiggle_smooth()} (R data frame with two columns: genome position and
#' smoothed signal).
#' @param wiggleData A data frame of wiggle data with two columns: genome position and signal. No default.
#' @param chr A number representing the chromosome of \code{wiggleData}. No default.
#' @param genome A string representing the genome used for mapping. No default.
#' @param yMax Optional number to be used as the max Y scale value in the plots.
#' Particularly useful to plot two chromosomes on the same Y scale. No default.
#' @param color Optional R color. Defaults to \code{grey50}.
#' @param protein A string representing the ChIPped protein. No default.
#' @param legendXcoord A number representing the X coordinate to locate legend.
#' Defaults to minimum X (left-aligned).
#' @param legendYcoord A number representing the Y coordinate to locate legend.
#' Defaults to maximum Y (top-aligned).
#' @param legendAnnotation Optional string to be used as the legend. Defaults to name of object
#' passed to the function.
#' @param onScreen Boolean indicating plots should be returned to the screen (\code{onScreen = TRUE})
#' or written to .png files (\code{onScreen = FALSE}). Defaults to \code{TRUE}.
#' @return A line plot, either on screen or as a png file (in the working directory).
#' @examples
#' \dontrun{
#' wiggle_plot(WT[[1]], 1, genome = 'SK1', protein = 'Red1')
#' 
#' wiggle_plot(WT_chr3, 3, genome = 'SK1', yMax = 5, color = 'red', protein = 'Red1', onScreen = TRUE)
#' 
#' wiggle_plot(chrXVI, 16, genome = 'S288C', yMax = 5, color = 'black',
#'             protein = 'Rec8-HA', legendXcoord = 600, onScreen = FALSE)
#' }
#' @export

wiggle_plot <- function(wiggleData, chr, genome, yMax, color = 'grey50', protein,
                        legendXcoord = -10, legendYcoord = yMax,
                        legendAnnotation = deparse(substitute(wiggleData)), onScreen = TRUE) {
  
  message('\nPlotting...')
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using script 'data-raw/data_internal.R'; stored in 'R/sysdata.rda'
  #----------------------------------------------------------------------------#
  # Load the data:
  if (genome == 'SK1') {
    Cen <- SK1cen
  } else {
    Cen <- S288Ccen
  }
  
  # Plot(s)  
  if (!onScreen) {
    png(filename = paste0(deparse(substitute(wiggleData)), "_chr", deparse(substitute(chr)),
                          ".png"), width = 1000, height = 480, units = 'px')
  }
  par(mfrow = c(1, 1), mar = c(8, 11, 2, 2), mgp = c(4, 2, 0))
  xMax <- ceiling(max(wiggleData[, 1], na.rm = T) / 1000)
  if (missing(yMax)) yMax <- ceiling(max(wiggleData[, 2],  na.rm = T))
  
  plot(wiggleData[, 1]/1000, wiggleData[, 2], type = 'l',
       lwd = 3, xaxt = 'n', yaxt = 'n',
       xlim = c(0, xMax),
       ylim = c(-2, yMax),
       xlab = paste0('Position on ', Cen[chr, 'Chromosome'], ' (Kb)'),
       ylab = paste0(protein, '\nChIP/Input'),
       main = paste0('Mapped to ', genome, ' genome'), col = color,
       bty = "n")
  axis(1, at = c(0, xMax))
  axis(2, at = c(0, yMax))
  
  points(Cen[chr, 4]/1000, -2.5, pch = 19)
  mtext(c('Cen'), 1, at = Cen[chr, 4]/1000, padj = 1)

  legend(legendXcoord, legendYcoord, legendAnnotation,
         bty = 'n', text.col = color)
  
  if (!onScreen) {
    dev.off()
  }
}