#' Average chromosome signal (coverage) dot plot
#'
#' This function allows you to plot average chromosome signal for one or two selected chromosomes.
#' It takes as input the output of 'chr_coverage' (16x2 R data frame).
#' @param coverageDataA A 16x2 data frame of coverage: chromosome and average signal. No default.
#' @param coverageDataB Optional 16x2 data frame of coverage: chromosome and average signal. No default.
#' @param protein A string representing the ChIPped protein. No default.
#' @param genome A string representing the genome used for mapping. Defaults to 'SK1'.
#' @param meanNorm Boolean indicating whether average coverage should be plotted as is
#' (meanNorm = FALSE) or normalized to genome-wide averages (meanNorm = TRUE). Defaults to TRUE.
#' @param yMax Optional number to be used as the max Y scale value in the plots. No default.
#' @param onscreen Boolean indicating plots should be returned to the screen (onScreen = TRUE)
#' or written to .png files (onScreen = FALSE). Defaults to TRUE.
#' @param colorA Optional R color for sample A. Defaults to 'grey50'.
#' @param colorB Optional R color for sample B. Defaults to 'green'.
#' @return A dot plot of one or two samples, either on screen or as a png file (in the working directory).
#' @examples
#' plot_chr_coverage(WT, rec8, protein = 'Red1', genome = 'SK1', meanNorm = TRUE, onScreen = TRUE, colorB = 'red')
#' 
#' plot_chr_coverage(WT, dot1, protein = 'Hop1', genome = 'S288C', meanNorm = FALSE, onScreen = FALSE)
#' @export

plot_chr_coverage <- function(coverageDataA, coverageDataB, protein,
                              genome = 'SK1', meanNorm = TRUE, yMax,
                              onScreen = TRUE,
                              colorA = 'grey50', colorB = 'green') {
  ptm <- proc.time()

  if (meanNorm) {
    ### Calculate coverage/chr relative to whole-genome average coverage
    # (i.e., normalize to whole-genome average)
    
    cat('Normalizing coverage/chr to genome average for',
        deparse(substitute(coverageDataA)), '\n')
    result <- as.data.frame(matrix(data = NA, nrow = 16, ncol = 2))
    result[, 1]  <-  coverageDataA[, 1]
    result[, 2]  <-  coverageDataA[, 2] / mean(coverageDataA[, 2])
    
    if (!missing(coverageDataB)) {
      cat('Normalizing coverage/chr to genome average for',
          deparse(substitute(coverageDataB)), '\n')
      result <- as.data.frame(matrix(data = NA, nrow = 16, ncol = 2))
      result[, 1]  <-  coverageDataB[, 1]
      result[, 2]  <-  coverageDataB[, 2] / mean(coverageDataB[, 2])
    }
  }
  
  ### Plot
  ##############################################################################
  # Information based on Keeney lab genome sequence and annotation
  # Not all this info is needed here; left it for future reference
  SK1cen <- data.frame("Chromosome" = c("chr01","chr02","chr03","chr04","chr05",
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
  # Information based on Wiggle Plots on Keeney lab genome sequence
  SK1rDNA <- data.frame("Name" = c("wiggle", "wiggle100"),
                        "Start" = c(433348, 433248), "End" = c(452165, 452265))
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
  data("S288Ccen", envir = parent.env(environment()))
  ##############################################################################
  
  cat('Plotting... \n')
  if (genome == 'SK1') {
    Cen <- SK1cen
  } else {
    Cen <- S288Ccen
  }
  
  lengths <- Cen[, c('Chromosome', 'LenChr')]
  ordered <- order(lengths[, 2])
  
  # Plot(s)  
  if (!onScreen) {
    png(filename = "Rplot.png", width = 1000, height = 480, unit = 'px')
  }
  
  
  if (!missing(coverageDataB)) {
    # Get highest of the maxima of the two strains to set y axis maximum
    if(missing(yMax)) {
      yMax_A <- ceiling(max(coverageDataA[, 2]))
      yMax_B <- ceiling(max(coverageDataB[, 2]))
      yMax <- tail(sort(c(yMax_A, yMax_B)), 1)
    }

  } else {
    if(missing(yMax)) yMax <- ceiling(max(coverageDataA[, 2]))
  }
  
  # Load package to use point transparency
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("R package 'scales' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  library(scales)
  
  par(mfrow = c(1, 1), mar = c(8, 12, 4, 2), mgp = c(6, 2, 0))
  plot(lengths[ordered, 2]/1000, coverageDataA[ordered, 2],
       xaxt = "n", yaxt = "n", xlim = c(0, 1500), ylim = c(-0.5, yMax),
       xlab = "Chromosome size (kb)", ylab = paste0(protein, '\nChIP/Input'),
       main = paste0('Mapped to ', genome, ' genome'),
       col = alpha(colorA, 0.7), pch = 19,
       cex = 3, cex.main = 2, cex.axis = 2, cex.lab = 2, bty = "n")
  axis(side = 1, at = c(0, 500, 1000, 1500), lwd = 4, cex.axis = 2, cex.lab = 2)
  axis(side = 2, at = c(0, 1, yMax), lwd = 4, cex.axis = 2, cex.lab = 2, las = 2)
  abline(h = 1, lty = 3, lwd = 2)
  
  if (!missing(coverageDataB)) {
    points(lengths[ordered, 2]/1000, coverageDataB[ordered, 2],
           col = alpha(colorB, 0.7), pch = 19, cex = 3)
  }
  
  if (missing(coverageDataB)) {
    legend(600, yMax, deparse(substitute(coverageDataA)), pch = 19,
           bty = 'n', pt.cex = 2, cex = 1.5,
           col = colorA, text.col = colorA)
  } else {
    legend(600, yMax,
           c(deparse(substitute(coverageDataA)), deparse(substitute(coverageDataB))),
           pch = 19, bty = 'n', pt.cex = 2, cex = 1.5,
           col = c(colorA, colorB), text.col = c(colorA, colorB)) 
  }
  
  if (!onScreen) dev.off()
  
  cat('... Completed in ', round((proc.time()[3] - ptm[3]), 2), ' sec.')
}