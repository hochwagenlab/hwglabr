#' Load wiggle data
#'
#' This function allows you to load tab-separated wiggle data.
#' Adapted from function written by Tovah Markowitz (original function name: 'readAll.tab').
#' @param fileLocation A string with the path to the folder containing the wiggle data. No default.
#' @param useReadr Boolean indicating whether to use the much faster 'read_tsv' from
#' Hadley Wickham's 'readr' package instead of base R's 'read.table'. Defaults to TRUE.
#' @param progressBar Boolean indicating whether to display a progress bar (using R package
#'  'pbapply'). Defaults to TRUE.
#' @param localCopy Boolean indicating whether to create local copy of target files before reading
#' data. If 'TRUE' a local folder is automatically created and deleted after use. Use this argument
#' to avoid reading files directly from shared locations (namely LabShare). Defaults to FALSE.
#' @return An R list of 16 data frames, one for each chromosome.
#' @examples
#' readall_tab("/Path/to/wiggles/folder/Red1_WT_reps_SacCer3_2mis_MACS_wiggle_norm/")
#' readall_tab("/Path/to/wiggles/", useReadr = TRUE, progressBar = TRUE)
#' readall_tab("/Path/to/wiggles/", useReadr = T, progressBar = T, localCopy = T)
#' @export

readall_tab <- function(fileLocation, useReadr = TRUE,
                        progressBar = TRUE, localCopy = FALSE) {
  ptm <- proc.time()
  if (localCopy) {
    cat('Copying data files to local folder "/temp"\n...')
    # Check if a directory 'temp' already exists
    if (file.exists('temp')) {
      stop('A folder called "temp" already exists in the current working directory.\n',
           'Please remove it and repeat function call.', call. = FALSE)
    }
    # Create temporary directory in current working directory
    # and make it the destination
    dir.create('temp')
    # Copy the files to the new temporary directory
    file.copy(fileLocation, 'temp', recursive = TRUE)
    # Update fileLocation to be the local directory
    fileLocation <- paste0('temp/', list.files('temp'))
  }
  
  cat('\nReading data\n')
  
  filenames <- list.files(fileLocation, full = T)
  if (length(filenames) == 17) {
    filenames <- filenames[2:17]
  }
  if (useReadr) {
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("R package 'readr' required. Please install it:\n", 
           "install.packages('readr')", call. = FALSE)
    }
    if (progressBar) {
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("R package 'pbapply' required. Please install it:\n", 
             "install.packages('pbapply')", call. = FALSE)
      }
      alldata <- pbapply::pblapply(filenames, readr::read_tsv, skip = 2,
                                   col_names = FALSE, progress = FALSE)
    }
    else {
      alldata <- lapply(filenames, readr::read_tsv, skip = 2,
                        col_names = FALSE, progress = FALSE)
    }
  }
  else {
    if (progressBar) {
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("R package 'pbapply' required. Please install it:\n", 
             "install.packages('pbapply')", call. = FALSE)
      }
      alldata <- pbapply::pblapply(filenames, read.table, skip = 2,
                                   sep = "\t")
    }
    else {
      alldata <- lapply(filenames, read.table, skip = 2, 
                        sep = "\t")
    }
  }
  
  if (localCopy) {
    cat('Deleting local copy\n...\n')
    # Delete temporary local directory
    unlink('temp', recursive = TRUE)
  }
  
  # Check which reference genome was used to map seq. data
  names(alldata) = filenames
  check_S288C <- any(grep("chrI.", names(alldata), fixed = TRUE))
  check_SK1 <- any(grep("chr01.", names(alldata), fixed = TRUE))
  if (check_S288C) {
    cat("Detected ref. genome - S288C\n(Chrs numbered using roman numerals)")
  }
  else if (check_SK1) {
    cat("\nDetected ref. genome - SK1\n(Chrs numbered using arabic numerals)")
  }
  else stop("Did not recognize reference genome.")
  
  cat("\n...\nCompleted in ", round((proc.time()[3] - ptm[3]), 1), " sec.")
  return(alldata)
}