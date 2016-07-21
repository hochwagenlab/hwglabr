#' Remove centromeric regions from wiggle data
#'
#' This function allows you to remove the region containing the centromere in each
#' chromosome. The original motivation was to analyse Red1 ChIP-seq data in rec8∆
#' mutants (show a strong enrichment around the centromeres).
#' @param wiggleData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.
#' @param regionSize Number indicating the size (in bp) of the region to remove
#' (centered on the centromere of each chromosome). Defaults to 25000 bp.
#' @return A similar list of the 16 chr wiggle data after removal of the rows for
#' the centromeric regions (both positions and signal) in each chromosome data frame.
#' @examples
#' remove_centromeres(rec8)
#' 
#' remove_centromeres(rec8, regionSize = 30000)
#' @export

remove_centromeres <- function(wiggleData, regionSize = 25000) {
  # Check reference genome 
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  check_S288C <- any(grep('chrI.', names(wiggleData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(wiggleData), fixed = TRUE))
  
  if (check_S288C) {
    cat('Ref. genome - S288C\n(Chrs numbered using roman numerals)\n')
    chrom <- chrom_S288C    
  } else if (check_SK1) {
    cat('Ref. genome - SK1\n(Chrs numbered using arabic numerals)\n')
    chrom <- chrom_SK1
  } else stop('Did not recognize reference genome.')
  
  cat('\nGetting centromeric locations...')
  #----------------------------------------------------------------------------#
  # All data loaded below is internal to the package
  # Generated using script 'data-raw/data_internal.R'; stored in 'R/sysdata.rda'
  #----------------------------------------------------------------------------#
  # Load the centromere data:
  load('R/sysdata.rda')
  
  if (check_SK1) {
    Cen <- SK1cen
  } else if(check_S288C) {
    Cen <- S288Ccen
  }
  
  cat('\nRemoving centromeric regions... \n')
  # Remove centromeric region data
  for (i in 1:length(wiggleData)) {
    chrNum <- paste0('chr', chrom[i], '.')
    
    # Get corresponding indices in the wiggle and centromere data
    wiggle_index <- grep(chrNum, names(wiggleData), fixed = TRUE)

    # Add '.' to chr # in Cen table to make them unique and then grep
    Cen$Chromosome <- paste0(Cen$Chromosome, ".")
    cen_index <- grep(chrNum, Cen$Chromosome, fixed = TRUE)
    
    ### Calculate positions centered on centromeres
    up_pos <- Cen[cen_index, 'Mid'] - (regionSize / 2)
    dwn_pos <- Cen[cen_index, 'Mid'] + (regionSize / 2)
    
    # Remove corresponding rows
    keep <- which(wiggleData[[wiggle_index]][, 1] <= up_pos | wiggleData[[wiggle_index]][, 1] >= dwn_pos)
    wiggleData[[wiggle_index]] <- wiggleData[[wiggle_index]][keep, ]
  }
  return (wiggleData)
}