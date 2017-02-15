#' Order chromosome column by chromosome length
#'
#' Given a data frame with 16 rows and a column with chromosome number, this function
#' returns the data frame with rows ordered by chromosome length.
#' @param inputData Data frame with chromosome column and 16 rows (for example the
#' output of \code{\link{chr_coverage}}).
#' No default.
#' @param chrColumn String name of chromosome column to order by. No default.
#' @param decreasing Logical indicating whether to order in decreasing order, from
#' longest to shortest chromosome. Defaults to \code{FALSE}.
#' @return Input data frame with rows ordered by chromosome length.
#' @examples
#' \dontrun{
#' ordered_df <- order_chromosomes(chrCov_df, 'seqname')
#' 
#' ordered_df <- order_chromosomes(chrCov_df, 'chr', decreasing=TRUE)
#' }
#' @export

order_chromosomes <- function(inputData, chrColumn, decreasing=FALSE) {
  # Make sure the input is a data frame with 16 rows
  if (!is.data.frame(inputData) | nrow(inputData) != 16) {
    stop("Wrong input data - not an R data frame wit 16 rows.", call. = FALSE)
  }
  
  # Check reference genome and get order
  chrom_S288C <- c('chrI', 'chrVI', 'chrIII', 'chrIX', 'chrVIII', 'chrV',
                   'chrXI', 'chrX', 'chrXIV', 'chrII', 'chrXIII', 'chrXVI',
                   'chrXII', 'chrVII', 'chrXV', 'chrIV')
  chrom_SK1 <- c('chr01', 'chr06', 'chr03', 'chr09', 'chr08', 'chr05', 'chr11',
                 'chr10', 'chr14', 'chr02', 'chr13', 'chr16', 'chr12', 'chr07',
                 'chr15', 'chr04')
  
  check_S288C <- any(grep('chrI', inputData[, chrColumn], fixed = TRUE))
  check_SK1 <- any(grep('chr01', inputData[, chrColumn], fixed = TRUE))
  
  if (check_S288C) {
    message('Ref. genome - S288C\n(Chrs numbered using roman numerals)')
    chrom <- chrom_S288C    
  } else if (check_SK1) {
    message('Ref. genome - SK1\n(Chrs numbered using arabic numerals)')
    chrom <- chrom_SK1
  } else stop('Did not recognize reference genome.\n',
              'Please make sure it is in one of the following two formats:\n',
              '"chrI" or "chr01"', call. = FALSE)
  
  if (decreasing) chrom <- rev(chrom)
  
  #inputData[[chrColumn]] <- factor(inputData[[chrColumn]], levels=chrom)
  return(inputData[match(chrom, inputData[[chrColumn]]), ])
}
