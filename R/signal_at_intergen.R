#' Signal at intergenic regions between convergent, divergent or tandem genes genome-wide
#'
#' This function allows you to pull out the ChIP signal over a window of positions of the selected size
#' (defaulting to \code{regionSize = 1000}) centered on midpoints of intergenic regions between genes of
#' the selected orientation: convergent, divergent or tandem.
#' It takes as input either the wiggle data as list of 16 chromosome (output of \code{\link{readall_tab}})
#' or the complete genome in one data frame (for example loaded from .bed files).\cr
#' 
#' @param inputData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}) or a data frame
#' (in which case you must set \code{inputDataFrame = TRUE}). No default.
#' @param inputDataFrame Boolean indicating whether input data is a data frame, as opposed to the
#' standard wiggle data in a list of 16 chromosomes. Defaults to \code{FALSE}.
#' @param orientation A string indicating the type of intergenic region to analyze, according to the
#' orientation of the flanking genes. Accepts one of:
#' \enumerate{
#'   \item \code{conv} Convergent genes
#'   \item \code{div} Divergent genes
#'   \item \code{tandem} Tandem genes. All regions in output will be in the same 5' -> 3' orientation
#'   (genes on Crick strand are reversed).
#' }
#' Defaults to \code{conv}.
#' @param regionSize Number indicating the size (in bp) of the region to calculate.
#' Defaults to 1000 bp (midpoint +/- 500 bp).
#' @param includeOverlapping Boolean indicating whether intergenic regions between overlapping genes should
#' be included in the analysis. If \code{includeOverlapping = FALSE}, intergenic regions between overlapping
#' genes will not be included in the analysis. Otherwise all regions are included. Defaults to \code{TRUE}.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current working
#' directory). If \code{saveFile = FALSE}, output is returned to screen or an R object (if assigned).
#' Defaults to \code{FALSE}.
#' @return A local data frame (dplyr data frame) with four columns:
#' \enumerate{
#'   \item \code{chr} Chromosome number
#'   \item \code{position} Nucleotide coordinate (relative to midpoint of intergenic region)
#'   \item \code{signal} ChIP signal
#'   \item \code{dist_apart} Distance between the two genes (negative value for overlapping genes)
#' }
#' @examples
#' signal_at_intergen(WT)
#' 
#' signal_at_intergen(WT, orientation = 'div', inputDataFrame = TRUE, regionSize = 1500,
#'                    includeOverlapping = FALSE, saveFile = TRUE)
#' @export

signal_at_intergen <- function(inputData, inputDataFrame = FALSE, orientation = 'conv',
                               regionSize = 1000, includeOverlapping = TRUE,
                               saveFile = FALSE) {
  ptm  <- proc.time()
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using 'data-raw/data_internal.R'; stored in 'R/sysdata.rda'
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
  
  # Check reference genome and load respective chromosome number vector
  if (check_S288C) {
    cat('Detected ref. genome - S288C\n')
    chrom <- chrom_S288C
  } else if (check_SK1) {
    print('Detected ref. genome - SK1')
    chrom <- chrom_SK1
  } else stop("Did not recognize reference genome.
              Check that chromosome numbers are in the usual format, e.g. 'chrI' or 'chr01'.")
  
  # Load intergenic region coordinate data
  if (check_S288C) {
    if (orientation == 'conv') {
      cat('Preparing convergent gene region info...\n')
      #data("S288C_conv_midpoint_dist")
      intergen <- S288C_conv_midpoint_dist
    } else if (orientation == 'div') {
      cat('Preparing divergent gene region info...\n')
      #data("S288C_div_midpoint_dist")
      intergen <- S288C_div_midpoint_dist
    } else if (orientation == 'tandem') {
      cat('Preparing tandem gene region info...\n')
      #data("S288C_tand_midpoint_dist")
      intergen <- S288C_tand_midpoint_dist
    } else stop("Did not recognize gene orientation.
                You should use one of the three strings 'conv', 'div' or 'tandem'.")
  } else if (check_SK1) {
    if (orientation == 'conv') {
      cat('Preparing convergent gene region info...\n')
      #data("SK1_conv_midpoint_dist")
      intergen <- SK1_conv_midpoint_dist
    } else if (orientation == 'div') {
      cat('Preparing divergent gene region info...\n')
      #data("SK1_div_midpoint_dist")
      intergen <- SK1_div_midpoint_dist
    } else if (orientation == 'tandem') {
      cat('Preparing tandem gene region info...\n')
      #data("SK1_tand_midpoint_dist")
      intergen <- SK1_tand_midpoint_dist
    } else stop("Did not recognize gene orientation.
                You should use one of the three strings 'conv', 'div' or 'tandem'.")
  }
  
  # Drop regions between overlapping genes?
  if (!includeOverlapping) {
    cat('Dropping overlapping genes...\n')
    intergen <- intergen[intergen$dist_apart > 0, ]
  }
  
  # Initialize object to collect final data
  allData <- data.frame()
  
  # Create data frame with intergen region info (midpoints +/- regionSize in bp)
  # If looking at tandem genes we need strand info
  if (orientation == 'tandem') {
    intergenicPos <- data.frame(matrix(nrow = nrow(intergen), ncol = 6))
    colnames(intergenicPos) <- c('chr', 'mid', 'up', 'dwn', 'dist_apart', 'strand')
  } else {
    intergenicPos <- data.frame(matrix(nrow = nrow(intergen), ncol = 5))
    colnames(intergenicPos) <- c('chr', 'mid', 'up', 'dwn', 'dist_apart')
  }
  intergenicPos[, 'chr'] <- intergen[, 'chr']                           # chr
  intergenicPos[, 'mid'] <- intergen[, 'midpoint']                      # mid
  intergenicPos[, 'up'] <- intergen[, 'midpoint'] - (regionSize / 2)    # up
  intergenicPos[, 'dwn'] <- intergen[, 'midpoint'] + (regionSize / 2)   # dwn
  intergenicPos[, 'dist_apart'] <- intergen[, 'dist_apart']             # overlap
  # If looking at tandem genes we need strand info in order to orient them all
  # 5' -> 3' at the end
  if (orientation == 'tandem') {
    intergenicPos[, 'strand'] <- intergen[, 'strand']                   # strand
  }
  
  cat('Collecting signal...\n')
  cat('(Skip any regions whose coordinates are not found in input data)\n\n')

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
    intergenGeneCount <- 0
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
      
      # For tandem genes in Crick strand, reverse the order of the position values
      # to have all regions oriented 5' -> 3'
      if (orientation == 'tandem' &&
          as.data.frame(intergenicPosChr)[j, 'strand'] == '-') {
          position <- (-1) * position
      }
      
      # Collect the data on the current intergenic region
      finalChromData <- dplyr::bind_rows(finalChromData,
                                         data.frame(chrNum, position, signal,
                                                    intergenicPosChr[j, 'dist_apart']))
      # Increment non-skipped gene count
      intergenGeneCount <- intergenGeneCount + 1
    }

    # Print number of skipped regions
    cat(paste0('... ', intergenGeneCount,
               ' intergenic regions (skipped ', j - intergenGeneCount, ')\n'))

    colnames(finalChromData) <- c('chr', 'position', 'signal', 'dist_apart')
    # Trim out NAs
    finalChromData <- finalChromData[complete.cases(finalChromData), ]
    # Collect the data on the current chromosome
    allData <- dplyr::bind_rows(allData, finalChromData)
  }
  
  cat(paste0('\nCompleted in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.\n'))
  
  if(saveFile) {
    cat('Saving file...\n')
    if(check_S288C) {
      write.table(allData, paste0(deparse(substitute(inputData)), "_", orientation,
                                  "_S288C.txt"), sep = "\t", quote = FALSE)
      cat('Done!')
    } else {
      write.table(allData, paste0(deparse(substitute(inputData)), "_", orientation,
                                  "_SK1.txt"), sep = "\t", quote = FALSE)
      cat('Done!')
    }
  } else {
    cat('Done!')
    return(allData)
  }
}