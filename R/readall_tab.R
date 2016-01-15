#' Load wiggle data
#'
#' This function allows you to load tab-separated wiggle data.
#' Adapted from function written by Tovah Markowitz (original function name: 'readAll.tab').
#' @param fileLocation A string with the path to the folder containing the wiggle data. No default.
#' @param use_readr Boolean indicating whether to use the much faster 'read_tsv' from
#' Hadley Wickham's 'readr' package instead of base R's 'read.table'. Defaults to FALSE.
#' @param progress_bar Boolean indicating whether to display a progress bar (using R package
#'  'pbapply'). Defaults to FALSE.
#' @return An R list of 16 data frames, one for each chromosome.
#' @examples
#' readall_tab("/Path/to/wiggles/folder/Red1_WT_reps_SacCer3_2mis_MACS_wiggle_norm/")
#' #' readall_tab("/Path/to/wiggles/", use_readr = TRUE, progress_bar = TRUE)
#' @export

readall_tab <- function(fileLocation, use_readr = FALSE, progress_bar = FALSE) {
  ptm <- proc.time()
  filenames <- list.files(fileLocation, full = T)
  # Check if 'all chromosomes' file exists and skipt it
  if (length(filenames) == 17) {
    filenames <- filenames[2:17]
  }
  
  if (use_readr) {
    # Check for package
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("R package 'readr' required. Please install it:\n",
           "install.packages('readr')", call. = FALSE)
    }
    library(readr)
    if (progress_bar) {
      # Check for package
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("R package 'pbapply' required. Please install it:\n",
             "install.packages('pbapply')", call. = FALSE)
      }
      library(pbapply)
      alldata <- pblapply(filenames, read_tsv, skip = 2, col_names = FALSE, progress = FALSE)
    } else {
      alldata <- lapply(filenames, read_tsv, skip = 2, col_names = FALSE, progress = FALSE)
    }
    
  } else {
    if (progress_bar) {
      # Check for package
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("R package 'pbapply' required. Please install it:\n",
             "install.packages('pbapply')", call. = FALSE)
      }
      library(pbapply)
      alldata <- pblapply(filenames, read.table, skip = 2, sep = "\t")
    } else {
      alldata <- lapply(filenames, read.table, skip = 2, sep = "\t")
    }
  }
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