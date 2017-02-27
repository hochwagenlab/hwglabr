#' Average chromosome signal (coverage) dot plot
#'
#' This function allows you to plot average chromosome signal for one or two selected samples.
#' It takes as input the output of \code{\link{chr_coverage}} (16x2 R data frame).
#' @param coverageDataA A 16x2 data frame of coverage: chromosome and average signal. No default.
#' @param coverageDataB Optional 16x2 data frame of coverage: chromosome and average signal. No default.
#' @param genome A string representing the genome used for mapping. No default.
#' @param meanNorm Boolean indicating whether input data are normalized to genome-wide
#' averages (\code{meanNorm = TRUE}). If so, a line is added at \code{y = 1}.
#' Defaults to FALSE.
#' @param yMax Optional number to be used as the max Y scale value in the plots. No default.
#' @param onScreen Boolean indicating plots should be returned to the screen (\code{onScreen = TRUE})
#' or written to .pdf files (\code{onScreen = FALSE}). Defaults to \code{TRUE}.
#' @param fileName A string to name the output .pdf file in case (\code{onScreen = FALSE}).
#' No default.
#' @param colorA Optional R color for sample A. Defaults to \code{grey50}.
#' @param colorB Optional R color for sample B. Defaults to \code{green}.
#' @return A dot plot of one or two samples, either on screen or as a .pdf file (in
#' the working directory).
#' @examples
#' \dontrun{
#' chr_coverage_plot(WT, rec8, genome = 'SK1', onScreen = TRUE, colorB = 'red')
#' 
#' chr_coverage_plot(WT, dot1, genome = 'S288C', meanNorm = TRUE,
#'                   onScreen = FALSE, fileName='chr_coverage_WT_and_dot1.pdf')
#' }
#' @export

chr_coverage_plot <- function(coverageDataA, coverageDataB, genome,
                              meanNorm = FALSE, yMax, onScreen = TRUE, fileName,
                              colorA = 'grey50', colorB = 'green') {
  ptm <- proc.time()
  
  # Make sure the input is the 16x2 data frame returned by chr_cov
  if (!is.data.frame(coverageDataA)) {
    stop("Wrong input data - not an R data frame.\n",
         "Please run 'chr_coverage' on your data first. Example:\n",
         "WT_cov_dataframe <- chr_coverage(WT_wiggle)\n",
         "chr_coverage_plot(WT_cov_dataframe, onScreen = TRUE)", call. = FALSE)
  } else {
    if (nrow(coverageDataA) != 16) {
      stop("Wrong input data dimensions.\n",
           "Please run 'chr_coverage' on your data first. Example:\n",
           "WT_cov_dataframe <- chr_coverage(WT_wiggle)\n",
           "chr_coverage_plot(WT_cov_dataframe, onScreen = TRUE)", call. = FALSE)
    }
  }
  
  if (!missing(coverageDataB)) {
    if (!is.data.frame(coverageDataB)) {
      stop(deparse(substitute(coverageDataB)), " is of wrong format - not an R data frame.\n",
           "Please run 'chr_coverage' on your data first. Example:\n",
           "WT_cov_dataframe <- chr_coverage(WT_wiggle)\n",
           "chr_coverage_plot(WT_cov_dataframe, genome = 'SK1', onScreen = TRUE)",
           call. = FALSE)
    } else {
      if (nrow(coverageDataB) != 16) {
        stop(deparse(substitute(coverageDataB)), " has wrong dimensions.\n",
             "Please run 'chr_coverage' on your data first. Example:\n",
             "WT_cov_dataframe <- chr_coverage(WT_wiggle)\n",
             "chr_coverage_plot(WT_cov_dataframe, genome = 'SK1', onScreen = TRUE)",
             call. = FALSE)
      }
    }
  }
  
  # Not changing variable names here causes problems with deparse(substitute())
  # in plot legend
  if (!missing(coverageDataB)) {
    dataA <- coverageDataA
    dataB <- coverageDataB
  } else {
    dataA <- coverageDataA
  }
  
  ### Plot
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using 'data-raw/data_internal.R'; stored in 'R/sysdata.rda'
  #----------------------------------------------------------------------------#
  
  # Load the data:
  message('Plotting... \n')
  if (genome == 'SK1') {
    Cen <- SK1cen
  } else {
    Cen <- S288Ccen
  }
  
  lengths <- Cen[, c('Chromosome', 'LenChr')]
  ordered <- order(lengths[, 2])
  
  # Plot(s)
  if (!onScreen) {
    #pdf(file = paste0(fileName, '_sizeBias.pdf'), width = 1000, height = 480, units = 'px')
    pdf(file = paste0(fileName), width = 4, height = 4)
  }
  
  if (!meanNorm) {
    col_to_plot <- 2
  } else {
    col_to_plot <- 4
  }
  
  if (!missing(coverageDataB)) {
    # Get highest of the maxima of the two strains to set y axis maximum
    if(missing(yMax)) {
      yMax_A <- ceiling(max(dataA[, col_to_plot]))
      yMax_B <- ceiling(max(dataB[, col_to_plot]))
      yMax <- tail(sort(c(yMax_A, yMax_B)), 1)
    }
  } else {
    if(missing(yMax)) yMax <- ceiling(max(dataA[, col_to_plot]))
  }
  
  if (!meanNorm) {
    ylab = 'ChIP-seq signal/Input'
  } else {
    ylab = 'ChIP-seq signal/Input\n(genome-wide average normalized)'
  }
  
  par(mar = c(5, 5, 4, 2))
  plot(lengths[ordered, 2]/1000, dataA[ordered, col_to_plot],
       xaxt = "n", yaxt = "n", xlim = c(0, 1550), ylim = c(-0.5, yMax),
       xlab = "Chromosome size (kb)", ylab = ylab,
       main = paste0('Mean signal per chromosome\nrelative to chromosome size'),
       col = colorA, pch = 19, cex = 1, cex.main = 1, cex.axis = 1, cex.lab = 1, bty = "n")
  axis(side = 1, at = c(0, 500, 1000, 1500), lwd = 1, cex.axis = 1, cex.lab = 1)
  axis(side = 2, at = c(0, 1, yMax), lwd = 1, cex.axis = 1, cex.lab = 1, las = 1)
  
  if (meanNorm) abline(h = 1, lty = 3, lwd = 1)
  
  if (!missing(coverageDataB)) {
    points(lengths[ordered, 2]/1000, dataB[ordered, col_to_plot],
           col = colorB, pch = 19, cex = 1)
  }
  
  if (!missing(coverageDataB)) {
    legend(600, yMax,
           c(deparse(substitute(coverageDataA)), deparse(substitute(coverageDataB))),
           pch = 19, bty = 'n', pt.cex = 1, cex = 1,
           col = c(colorA, colorB), text.col = c(colorA, colorB))
  }
  
  if (!onScreen) dev.off()
  
  message('... Completed in ', round((proc.time()[3] - ptm[3]), 2), ' sec.')
}