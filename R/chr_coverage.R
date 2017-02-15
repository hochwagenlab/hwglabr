#' Average signal per chromosome
#'
#' This function allows you to calculate the average signal (or coverage) per chromosome.
#' @param wiggleData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}).
#' No default.
#' @param meanNorm Boolean indicating whether to normalize by average genome-wide
#' signal. This adds two columns to output: mean-subtracted signal and mean-divided
#' signal. Defaults to \code{FALSE}.
#' @param orderChrs Boolean indicating whether to order rows by chromosome length.
#' Defaults to \code{FALSE}.
#' @param removeCen Boolean indicating whether to remove regions around centromeres
#' (using function \code{\link{remove_centromeres}}). Defaults to \code{FALSE}.
#' @param cenRegionSize Number indicating the size (in bp) of the region to remove
#' (centered on the centromere of each chromosome). Corresponds to argument \code{regionSize}
#' of \code{\link{remove_centromeres}}. Defaults to 50'000 bp.
#' @return A 16x2 R data frame with either two or four columns:
#' \enumerate{
#'   \item \code{chr} Chromosome number
#'   \item \code{mean_coverage} average signal
#'   \item \code{mean_sub_coverage} signal minus genome-wide average (chr_signal - genome_signal)
#'   \item \code{mean_div_coverage} signal divided by genome-wide average (chr_signal / genome_signal)
#' }
#' @examples
#' \dontrun{
#' chr_coverage(WT)
#' 
#' chr_coverage(dot1, removeCen = TRUE, cenRegionSize = 50000)
#' }
#' @export

chr_coverage <- function(wiggleData, meanNorm=FALSE, orderChrs=FALSE,
                         removeCen = FALSE, cenRegionSize = 50000) {
  # Make sure the input is a list with 16 elements
  if (!is.list(wiggleData) | length(wiggleData) != 16) {
    stop("Wrong input data - not a list with 16 elements.\n",
         "Please provide a list of 16 data frames, one for each of the 16 chromosomes",
         "(output of readall_tab).", call. = FALSE)
  }
  
  # Check reference genome 
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  check_S288C <- any(grep('chrI.', names(wiggleData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(wiggleData), fixed = TRUE))
  
  if (check_S288C) {
    message('Ref. genome - S288C\n(Chrs numbered using roman numerals)')
    chrom <- chrom_S288C    
  } else if (check_SK1) {
    message('Ref. genome - SK1\n(Chrs numbered using arabic numerals)')
    chrom <- chrom_SK1
  } else stop('Did not recognize reference genome.')
  
  if (meanNorm) {
    coverageTable <- as.data.frame(matrix(0, nrow = length(wiggleData), ncol = 4))
    gen_avrg <- suppressMessages(hwglabr::genome_average(wiggleData))
  } else {
    coverageTable <- as.data.frame(matrix(0, nrow = length(wiggleData), ncol = 2))
  }
  
  if(removeCen){
    wiggleData <- hwglabr::remove_centromeres(wiggleData = wiggleData,
                                              regionSize = cenRegionSize)
  }
  
  for (i in 1:length(wiggleData)) {
    chrNum <- paste0('chr', chrom[i], '.')
    index <- grep(chrNum, names(wiggleData), fixed = TRUE)
    coverageTable[i, 1] <- paste0('chr', chrom[i])
    # mean() function does not work if dplyr was loaded when the wiggle data was
    # loaded (because of dplyr's non standard evaluation: data is in tbl_df class)
    # Workaround: convert to data frame before passing to mean()
    coverageTable[i, 2] <- mean(as.data.frame(wiggleData[[index]])[, 2], na.rm=TRUE)
    if (meanNorm) {
      coverageTable[i, 3] <- coverageTable[i, 2] - gen_avrg
      coverageTable[i, 4] <- coverageTable[i, 2] / gen_avrg
    }
  }
  
  if (meanNorm) {
    colnames(coverageTable) <- c("chr", "mean_coverage",
                                 "mean_sub_coverage", "mean_div_coverage")
  } else {
    colnames(coverageTable) <- c("chr", "mean_coverage")
  }
  
  if (orderChrs) coverageTable <- hwglabr::order_chromosomes(coverageTable, 'chr')
  
  return (coverageTable)
}
