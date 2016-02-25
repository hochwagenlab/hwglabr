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
#' @param onscreen Boolean indicating plots should be returned to the screen (\code{onScreen = TRUE})
#' or written to .png files (\code{onScreen = FALSE}). Defaults to \code{TRUE}.
#' @return A line plot, either on screen or as a png file (in the working directory).
#' @examples
#' wiggle_plot(WT[[1]], 1, genome = 'SK1', protein = 'Red1')
#' 
#' wiggle_plot(WT_chr3, 3, genome = 'SK1', yMax = 5, color = 'red', protein = 'Red1', onScreen = TRUE)
#' 
#' wiggle_plot(chrXVI, 16, genome = 'S288C', yMax = 5, color = 'black',
#' protein = 'Rec8-HA', legendXcoord = 600, onScreen = FALSE)
#' @export

wiggle_plot <- function(wiggleData, chr, genome, yMax, color = 'grey50', protein,
                        legendXcoord = -10, legendYcoord = yMax, onScreen = TRUE) {
  
  ##############################################################################
  # Information based on Keeney lab genome sequence and annotation
  SK1cen<- data.frame("Chromosome" = c("chr01","chr02","chr03","chr04","chr05",
                                       "chr06","chr07","chr08","chr09","chr10",
                                       "chr11","chr12","chr13","chr14","chr15",
                                       "chr16"),
                      "Start" = c(137832, 226711, 128699, 463204, 157003, 162815,
                                  505440, 95031, 346215, 415648, 452723, 137738,
                                  249103, 616840, 307236, 553355),
                      "End" = c(137948, 226826, 128779, 463321, 157119, 162931,
                                505558, 95147, 346330, 415764, 452838, 137855,
                                249221, 616956, 307353, 553467),
                      "Mid" = c(137890, 226768, 128739, 463262, 157061, 162873,
                                505499, 95089, 346272, 415706, 452780, 137796,
                                249162, 616898, 307294, 553411),
                      "LenChr" = c(203893, 794508, 342718, 1490682, 602514,
                                   284456, 1067526, 544538, 435585, 719294,
                                   687260, 1008248, 908607, 812465, 1054033,
                                   921188))

  ##############################################################################
  ##############################################################################
  # Generate the same information for S288C genome
  # The data is internal to the package and was generated as follows:
  #
  # path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/'
  # S288C_gff <- read.table(paste0(path, 'saccharomyces_cerevisiae_R64-1-1_20110208.gff'),
  #                         fill = TRUE, stringsAsFactors = FALSE)
  # S288C_gff <- S288C_gff[1:16406, ]
  # S288Ccen <- S288C_gff[S288C_gff[, 3] == 'centromere', c(1, 4:5)]
  # names(S288Ccen) <- c('Chromosome', 'Start', 'End')
  # S288Ccen$Mid <- floor(S288Ccen$Start +  (S288Ccen$End - S288Ccen$Start) / 2)
  # S288Ccen$LenChr <- S288C_gff[S288C_gff[, 3] == 'chromosome', 5][1:16]
  # setwd('/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/Scripts/Rpackages/hwglabr')
  # devtools::use_data(S288Ccen, internal = FALSE)
  
  # Load the data:
  data("S288Ccen")
  ##############################################################################
  
  cat('Plotting... \n')
  if (genome == 'SK1') {
    Cen <- SK1cen
  } else {
    Cen <- S288Ccen
  }
  
  # Plot(s)  
  if (!onScreen) {
    png(filename = paste0(deparse(substitute(wiggleData)), "_chr", deparse(substitute(chr)),
                          ".png"), width = 1000, height = 480, unit = 'px')
  }
  par(mfrow = c(1, 1), mar = c(8, 11, 2, 2), mgp = c(6, 2, 0))
  xMax <- ceiling(max(wiggleData[, 1]) / 1000)
  if (missing(yMax)) yMax <- ceiling(max(wiggleData[, 2]))
  
  plot(wiggleData[, 1]/1000, wiggleData[, 2], type = 'l',
       lwd = 3, xaxt = 'n', yaxt = 'n',
       xlim = c(0, xMax),
       ylim = c(-2, yMax),
       xlab = paste0('Chr', chr, ' position (Kb)'),
       ylab = paste0(protein, '\nChIP/Input'),
       main = paste0('Mapped to ', genome, ' genome'), col = color,
       cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, bty = "n")
  axis(1, at = c(0, xMax), lwd = 4, cex.axis = 2, cex = 2.5)
  axis(2, at = c(0, yMax), lwd = 4, las = 2, cex.axis = 2, cex = 2.5)
  
  points(Cen[chr, 4]/1000, -2.5, pch = 19, cex = 3.0)
  mtext(c('Cen'), 1, at = Cen[chr, 4]/1000, cex = 1.5, padj = 1)

  legend(legendXcoord, legendYcoord, deparse(substitute(wiggleData)),
         bty = 'n', cex = 3, text.col = color)
  
  if (!onScreen) {
    dev.off()
  }
}