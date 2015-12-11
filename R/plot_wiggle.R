#' Wiggle data line plot
#'
#' This function allows you to plot wiggle data for one or two selected chromosomes.
#' @param wiggleDataA A data frame of wiggle data with two columns: genome position and signal. No default.
#' @param wiggleDataB Optional data frame of wiggle data with two columns: genome position and signal. No default.
#' @param chrA A number representing the chromosome of 'wiggleDataA'. No default.
#' @param chrB Optional number representing the chromosome of 'wiggleDataB'.
#' Must be provided if 'wiggleDataB' is used. No default.
#' @param genome A string representing the genome used for mapping. Defaults to 'SK1'.
#' @param yMax Optional number to be used as the max Y scale value in the plots.
#' Particularly useful to plot two chromosomes on the same Y scale. No default.
#' @param color Optional R color. Defaults to 'grey50'.
#' @param protein A string representing the ChIPped protein. No default.
#' @param onscreen Boolean indicating plots should be returned to the screen (onScreen = TRUE)
#' or written to .png files (onScreen = FALSE). Defaults to TRUE.
#' @return One or two line plots, either on screen or as png files (in the working directory).
#' @examples
#' plot_wiggle(WT_chr3, WT_chr5, 3, 5, genome = 'SK1', yMax = 5, color = 'red', protein = 'Red1', onScreen = TRUE)
#' 
#' plot_wiggle(chr1, chr10, 1, 10, genome = 'S288C', yMax = 5, color = 'black', protein = 'Rec8-HA', onScreen = FALSE)

plot_wiggle <- function(wiggleDataA, wiggleDataB, chrA, chrB, genome = 'SK1',
                        yMax, color = 'grey50', protein, onScreen = TRUE) {
  ptm <- proc.time()
  
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
  # devtools::use_data(S288Ccen, internal = TRUE)
  
  # Load the internal data:
  data(sysdata, envir = environment())
  
  ##############################################################################
  
  cat('Plotting... ')
  if (genome == 'SK1') {
    Cen <- SK1cen
  } else {
    Cen <- S288Ccen
  }
  
  # Plot(s)  
  if (!onScreen) {
    png(filename = "RplotA.png", width = 1000, height = 480, unit = 'px')
  }
  par(mfrow = c(1, 1), mar = c(8, 11, 2, 2), mgp = c(6, 2, 0))
  xMaxA <- ceiling(max(wiggleDataA[, 1]) / 1000)
  ifelse(missing(yMax), yMaxA <- ceiling(max(wiggleDataA[, 2])), yMaxA <- yMax)
  
  plot(wiggleDataA[, 1]/1000, wiggleDataA[, 2], type = 'l',
       lwd = 3, xaxt = 'n', yaxt = 'n',
       xlim = c(0, xMaxA),
       ylim = c(-2, yMaxA),
       xlab = paste0('Chr', chrA, ' position (Kb)'),
       ylab = paste0(protein, '\nChIP/Input'),
       main = paste0('Mapped to ', genome, ' genome'), col = color,
       cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, bty = "n")
  axis(1, at = c(0, xMaxA), lwd = 4, cex.axis = 2, cex = 2.5)
  axis(2, at = c(0, yMaxA), lwd = 4, las = 2, cex.axis = 2, cex = 2.5)
  
  points(Cen[chrA, 4]/1000, -2.5, pch = 19, cex = 3.0)
  mtext(c('Cen'), 1, at = Cen[chrA, 4]/1000, cex = 1.5, padj = 1)
  #abline(h = 0, lty = 3)
  legend(-10, yMaxA, deparse(substitute(wiggleDataA)),
         bty = 'n', cex = 3, text.col = color)
  
  if (!onScreen) {
    dev.off()
  }

  # Handling a second chromosome -----------------------------------------------
  if(!missing(wiggleDataB)) {
    if (!onScreen) {
      png(filename = "RplotB.png", width = 1000, height = 480, unit = 'px')
    }
    par(mfrow = c(1, 1), mar = c(8, 11, 2, 2), mgp = c(6, 2, 0))
    xMaxB <- ceiling(max(wiggleDataB[, 1]) / 1000)
    ifelse(missing(yMax), yMaxB <- ceiling(max(wiggleDataB[, 2])), yMaxB <- yMax)
    
    plot(wiggleDataB[, 1]/1000, wiggleDataB[, 2], type = 'l',
         lwd = 3, xaxt = 'n', yaxt = 'n',
         xlim = c(0, xMaxB),
         ylim = c(-2, yMaxB),
         xlab = paste0('Chr', chrB, ' position (Kb)'),
         ylab = paste0(protein, '\nChIP/Input'),
         main = paste0('Mapped to ', genome, ' genome'), col = color,
         cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, bty = "n")

    axis(1, at = c(0, xMaxB), lwd = 4, cex.axis = 2, cex = 2.5)
    axis(2, at = c(0, yMaxB), lwd = 4, las = 2, cex.axis = 2, cex = 2.5)
      
    points(Cen[chrB, 4]/1000, -2.5, pch = 19, cex = 3.0)
    mtext(c('Cen'), 1, at = Cen[chrB, 4]/1000, cex = 1.5, padj = 1)
    #abline(h = 0, lty = 3)
    legend(-10, yMaxB, deparse(substitute(wiggleDataB)),
           bty = 'n', cex = 3, text.col = color)
    
    if (!onScreen) {
      dev.off()
    }
  }

  cat('\n...\nCompleted in', round((proc.time()[3] - ptm[3]), 1), 'sec.', sep = ' ')
}