#' Signal between all convergent genes genomewide
#'
#' This function allows you to pull out the ChIP signal centered on midpoints of convergent genes.
#' It takes as input either the wiggle data as list of 16 chromosome (output of 'readall_tab()')
#' or complete genome in one data frame (for example loaded from .bed files).
#' @param inputData As a list of the 16 chr wiggle data (output of readall_tab). No default.
#' @param regionSize Number indicating the size (in bp) of the region to calculate.
#' Defaults to 1000 bp (+/- 500 bp).
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current working
#' directory). If 'saveFile = FALSE', output is returned to screen or an R object (if assigned).
#' Defaults to FALSE.
#' @param input_dataFrame Boolean indicating whether input data is a data frame. This is the case when you
#' have data loaded from a .bed format, for example, typically nucleosome signal, as opposed to the standard
#' wiggle data in a list of 16 chromosomes. Defaults to FALSE.
#' @return An R data frame with two columns: genome position and mean signal.
#' @examples
#' signal_at_conv(WT)
#' 
#' signal_at_conv(WT, region_size = 1500, saveFile = TRUE, input_dataFrame = FALSE)
#' @export

signal_at_conv <- function(inputData, region_size = 1000, saveFile = FALSE,
                           input_dataFrame = FALSE) {
  ptm  <- proc.time()
  
  ##############################################################################
  # Get the convergent gene region data
  # The data is internal to the package and was generated as follows:
  #
  #-------------------------------- S288C genome ------------------------------#
  # Import data file, originally generated using:
  # "/Volumes/LabShare/Luis/LabWork/Scripts/2015.10_H3k79me3_transcription.R"
  # path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/2015.05_dot1_rec8dot1/2015.10_H3K79me3_transcripts/'
  # S288C_conv <- read.table(paste0(path, 'conv_midpoints_rpkm.txt'),
  #                         header = TRUE, stringsAsFactors = FALSE)
  # Remove info about transcription:
  # S288C_conv <- S288C_conv[, 1:6]
  # setwd('/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/Scripts/Rpackages/hwglabr')
  # devtools::use_data(S288C_conv, internal = FALSE)
  #
  #--------------------------------- SK1 genome -------------------------------#
  # Import data file, originally generated using:
  # /Volumes/LabShare/Luis/LabWork/Scripts/2015.11_convMidpointsSK1.R
  # path <- '/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/GenomeSequences/SK1/'
  # SK1_conv <- read.table(paste0(path, 'conv_midpoints_SK1.txt'),
  #                        header = TRUE, stringsAsFactors = FALSE)
  # setwd('/Users/luis/Google_Drive_NYU/LabShare_Luis/LabWork/Scripts/Rpackages/hwglabr')
  # devtools::use_data(SK1_conv, internal = FALSE)
  #
  
  # Check reference genome and load appropriate convergent gene regions 
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  check_S288C <- any(grep('chrI.', names(inputData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(inputData), fixed = TRUE))
  
  if (check_S288C) {
    cat('Detected ref. genome - S288C\n')
    # Load the data:
    data("S288C_conv")
    conv <- S288C_conv
    chrom <- chrom_S288C
    
  } else if (check_SK1) {
    print('Detected ref. genome - SK1')
    # Load the data:
    data("SK1_conv")
    conv <- SK1_conv
    chrom <- chrom_SK1
  } else stop('Did not recognize reference genome.')
  
  cat('Preparing convergent gene region info...\n')
  
  # Initialize object to collect final data
  allData <- data.frame()
  
  # Create data frame with conv region info (midpoints +/- region_size in bp)
  intergenicPos <- data.frame(matrix(nrow = nrow(conv), ncol = 4))
  colnames(intergenicPos) <- c('chr', 'mid', 'up', 'dwn')
  intergenicPos[, 'chr'] <- conv[, 'chr']                            # chr
  intergenicPos[, 'mid'] <- conv[, 'midpoint']                       # mid
  intergenicPos[, 'up'] <- conv[, 'midpoint'] - (region_size / 2)    # up
  intergenicPos[, 'dwn'] <- conv[, 'midpoint'] + (region_size / 2)   # dwn
  
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  library(dplyr)
  
  cat('Collecting signal...\n')
  
  for(i in 1:length(chrom)) {
    chrNum <- paste0('chr', chrom[i])
    # Get ChIP data list item corresponding to chrom to analyze
    # Either list of wiggles or dataFrame from .bed ('input_dataFrame = F / T')
    if (input_dataFrame) {
      chromData <- inputData[inputData[, 1] == chrom[i], ]
    } else {
      # Get index of ChIP data list item corresponding to chrom to analyze
      # Add '.' to make it unique (otherwise e.g. 'chrI' matches 'chrII' too)
      listIndex <- grep(paste0(chrNum, '.'), names(inputData), fixed = TRUE)
      chromData <- inputData[[listIndex]]
      colnames(chromData) <- c('position', 'signal')
    }

    finalChromData <- data.frame()

    # Get subset of intergenic regions on chr i
    intergenicPosChr <- intergenicPos[intergenicPos[, 'chr'] == chrNum, ]
    
    for(j in 1:nrow(intergenicPosChr)) {
      if(j == 1) {
        cat(paste0('\n', chrNum, ':\n'))
      }
      
      if(j == nrow(intergenicPosChr)) {
        cat(paste0('\n', j, ' intergenic regions\n'))
      }
      # Skip if gene coordinates not in ChIPseq data
      #      if(!intergenicPosChr[j, 'up'] %in% chromData[, 'position'] |
      #         !intergenicPosChr[j, 'dwn'] %in% chromData[, 'position']) {
      if(!any(chromData[, 'position'] == intergenicPosChr[j, 'up']) |
         !any(chromData[, 'position'] == intergenicPosChr[j, 'dwn'])) {
        next
      }
      
      upPosition <- intergenicPosChr[j, 'up']
      dwnPosition <- intergenicPosChr[j, 'dwn']
      # position
      position <- chromData[chromData[, 'position'] >= upPosition &
                              chromData[, 'position'] <= dwnPosition, 'position']
      position <- position - intergenicPosChr[j, 'mid']
      # signal
      signal <- chromData[chromData[, 'position'] >= upPosition &
                            chromData[, 'position'] <= dwnPosition, 'signal']
      
      # collect the data
      finalChromData <- bind_rows(finalChromData, data.frame(chrNum, position, signal))
    }
    colnames(finalChromData) <- c('chr', 'position', 'signal')
    # Trim out NAs
    finalChromData <- finalChromData[complete.cases(finalChromData), ]
    allData <- bind_rows(allData, finalChromData)
  }
  if(saveFile) {
    cat(paste0('Saving file...\n'))
    if(check_S288C) {
      write.table(allData, paste0(deparse(substitute(inputData)),
                                  "_S288C_conv.txt"), sep = "\t", quote = FALSE)  
    } else {
      write.table(allData, paste0(deparse(substitute(inputData)),
                                  "_SK1_conv.txt"), sep = "\t", quote = FALSE)
    }
  } else {
    return(allData)
  }
  cat(paste0('Completed in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.\n'))
}