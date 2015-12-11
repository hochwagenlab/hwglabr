#' Load wiggle data
#'
#' This function allows you to load tab-separated wiggle data.
#' @param fileLocation The path to the folder containing the wiggle data. No default.
#' @return An R list of 16 data frames, one for each chromosome.
#' @examples
#' readall_tab(/Path/to/wiggle/files/folder)

# To load ChIP-seq tab-delimited wiggle files #
readall_tab <- function(fileLocation) {
  ptm <- proc.time()
  filenames <- list.files(fileLocation, full = T)
  # check if 'all chromosomes' file exists
  if (length(filenames) == 17) {
    filenames <- filenames[2:17]
  }  
  alldata <- lapply(filenames, read.table, skip = 2, sep = "\t")
  names(alldata) = filenames
    
  # Check reference genome
  check_S288C <- any(grep('chrI.', names(alldata), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(alldata), fixed = TRUE))
  
  if (check_S288C) {
    cat('Detected ref. genome - S288C\n(Chrs numbered using roman numerals)')
  } else if (check_SK1) {
    cat('Detected ref. genome - SK1\n(Chrs numbered using arabic numerals)')
  } else stop('Did not recognize reference genome.')
  
  cat('\n...\nCompleted in', round((proc.time()[3] - ptm[3]), 1), 'sec.', sep = ' ')
  return (alldata)
}