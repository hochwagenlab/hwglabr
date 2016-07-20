#' Average signal per chromosome
#'
#' This function allows you to calculate the average signal (or coverage) per chromosome.
#' @param wiggleData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.
#' @return A 16x2 R data frame with two columns: chromosome number and average signal.
#' @examples
#' chr_coverage(WT)
#' 
#' chr_coverage(rec8)
#' @export

chr_coverage <- function(wiggleData) {
  coverageTable <- as.data.frame(matrix(0, nrow = length(wiggleData), ncol = 2))
  
  # Check reference genome 
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  check_S288C <- any(grep('chrI.', names(wiggleData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(wiggleData), fixed = TRUE))
  
  if (check_S288C) {
    cat('Ref. genome - S288C\n(Chrs numbered using roman numerals)')
    chrom <- chrom_S288C    
  } else if (check_SK1) {
    cat('Ref. genome - SK1\n(Chrs numbered using arabic numerals)')
    chrom <- chrom_SK1
  } else stop('Did not recognize reference genome.')
  
  for (i in 1:length(wiggleData)) {
    chrNum <- paste0('chr', chrom[i], '.')
    index <- grep(chrNum, names(wiggleData), fixed = TRUE)
    coverageTable[i, 1] <- paste0('chr', chrom[i])
    # mean() function does not work if dplyr was loaded when the wiggle data was loaded
    # (because of dplyr's non standard evaluation: data is in tbl_df class)
    # Workaround: write out division of sum of values by their number in case
    # mean() is not working
    coverageTable[i, 2] <- mean(wiggleData[[index]][, 2], na.rm = T)
    if(is.na(mean(wiggleData[[index]][, 2]))){
      coverageTable[i, 2] <- sum(wiggleData[[index]][, 2]) / nrow(wiggleData[[index]][, 2]) 
    }
  }
  colnames(coverageTable) <- c("chr", "mean_coverage")
  return (coverageTable)
}