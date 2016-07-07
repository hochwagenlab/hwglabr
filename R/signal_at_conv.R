#' Signal between all convergent genes genome-wide
#'
#' \strong{Deprecated! Use \code{\link{signal_at_intergen}} instead.} \cr
#' This function allows you to pull out the ChIP signal centered on midpoints of convergent genes.
#' It takes as input either the wiggle data as list of 16 chromosome (output of \code{\link{readall_tab}})
#' or complete genome in one data frame (for example loaded from .bed files).
#' @param inputData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}) or a data frame
#' (in which case you must set \code{inputDataFrame = TRUE}). No default.
#' @param regionSize Number indicating the size (in bp) of the region to calculate.
#' Defaults to 1000 bp (+/- 500 bp).
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current working
#' directory). If \code{saveFile = FALSE}, output is returned to screen or an R object (if assigned).
#' Defaults to \code{FALSE}.
#' @param inputDataFrame Boolean indicating whether input data is a data frame. This is the case when you
#' have data loaded from a .bed format, typically nucleosome signal, as opposed to the standard wiggle data
#' in a list of 16 chromosomes. Defaults to \code{FALSE}.
#' @return A local data frame (dplyr data frame) with three columns: chr (chromosome number), position
#' (genome coordinate relative to midpoint of intergenic region) and signal (ChIP signal).
#' @examples
#' signal_at_conv(WT)
#' 
#' signal_at_conv(WT, regionSize = 1500, saveFile = TRUE, inputDataFrame = FALSE)
#' @export

signal_at_conv <- function(inputData, regionSize = 1000, saveFile = FALSE,
                           inputDataFrame = FALSE) {
  ptm  <- proc.time()
  
  cat('Note: This function is deprecated!\n',
      'From version 0.2 on it is no longer maintained. Use signal_at_intergen() instead.')
  
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using 'data-raw/data_internal.R'; stored in 'data/sysdata.rda'
  #----------------------------------------------------------------------------#
  
  # Check reference genome and load appropriate convergent gene regions 
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  
  # Handle the case 'inputDataFrame = T'
  if (inputDataFrame) {
    check_S288C <- any(inputData[, 1] == 'chrI')
    check_SK1 <- any(inputData[, 1] == 'chr01')
  } else {
    check_S288C <- any(grep('chrI.', names(inputData), fixed = TRUE))
    check_SK1 <- any(grep('chr01.', names(inputData), fixed = TRUE))
  }
  
  # Load the data:
  data(sysdata, envir=environment())
  
  if (check_S288C) {
    cat('Detected ref. genome - S288C\n')
    conv <- S288C_conv
    chrom <- chrom_S288C
    
  } else if (check_SK1) {
    print('Detected ref. genome - SK1')
    conv <- SK1_conv
    chrom <- chrom_SK1
  } else stop("Did not recognize reference genome.
              Check that chromosome numbers are in the usual format, e.g. 'chrI' or 'chr01'.")
  
  cat('Preparing convergent gene region info...\n')
  
  # Initialize object to collect final data
  allData <- data.frame()
  
  # Create data frame with conv region info (midpoints +/- regionSize in bp)
  intergenicPos <- data.frame(matrix(nrow = nrow(conv), ncol = 4))
  colnames(intergenicPos) <- c('chr', 'mid', 'up', 'dwn')
  intergenicPos[, 'chr'] <- conv[, 'chr']                           # chr
  intergenicPos[, 'mid'] <- conv[, 'midpoint']                      # mid
  intergenicPos[, 'up'] <- conv[, 'midpoint'] - (regionSize / 2)    # up
  intergenicPos[, 'dwn'] <- conv[, 'midpoint'] + (regionSize / 2)   # dwn
  
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  cat('Collecting signal...\n')
  cat('(Skip any regions whose coordinates are not found in input data)\n')

  for(i in 1:length(chrom)) {
    chrNum <- paste0('chr', chrom[i])
    # Get ChIP data list item corresponding to chrom to analyze
    # Either list of wiggles or dataFrame ('inputDataFrame = F / T')
    if (inputDataFrame) {
      chromData <- inputData[inputData[, 1] == paste0('chr', chrom[i]), ]
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
    
    # Count intergenic regions skipped due to missing coordinates in wiggle data
    convGeneCount <- 0
    for(j in 1:nrow(intergenicPosChr)) {
      if(j == 1) {
        cat(paste0(chrNum, ':\n'))
      }
      
      # Skip if gene coordinates not in ChIPseq data
      # Comparison below will not work if dplyr was loaded when the wiggle data was loaded
      # (because of dplyr's non standard evaluation: data is in tbl_df class)
      # Workaround: convert to data.frame first
      #      if(!intergenicPosChr[j, 'up'] %in% chromData[, 'position'] |
      #         !intergenicPosChr[j, 'dwn'] %in% chromData[, 'position']) {
      if(!any(as.data.frame(chromData)[, 'position'] == as.data.frame(intergenicPosChr)[j, 'up']) |
         !any(as.data.frame(chromData)[, 'position'] == as.data.frame(intergenicPosChr)[j, 'dwn'])) {
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
      finalChromData <- dplyr::bind_rows(finalChromData, data.frame(chrNum, position, signal))
      # Increment non-skipped gene count
      convGeneCount <- convGeneCount + 1
    }

    # Print number of skipped regions
    cat(paste0('... ', convGeneCount,
               ' intergenic regions (skipped ', j - convGeneCount, ')\n'))

    colnames(finalChromData) <- c('chr', 'position', 'signal')
    # Trim out NAs
    finalChromData <- finalChromData[complete.cases(finalChromData), ]
    allData <- dplyr::bind_rows(allData, finalChromData)
  }
  
  cat(paste0('\nCompleted in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.\n'))
  
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
}