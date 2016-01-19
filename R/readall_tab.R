#' Load wiggle data
#'
#' This function allows you to load tab-separated wiggle data.
#' Adapted from function written by Tovah Markowitz (original function name: 'readAll.tab').
#' @param fileLocation A string with the path to the folder containing the wiggle data. No default.
#' @param use_readr Boolean indicating whether to use the much faster 'read_tsv' from
#' Hadley Wickham's 'readr' package instead of base R's 'read.table'. Defaults to FALSE.
#' @param progress_bar Boolean indicating whether to display a progress bar (using R package
#'  'pbapply'). Defaults to FALSE.
#' @param local_copy Boolean indicating whether to create local copy of target files before reading
#' data. If 'TRUE' a local folder is automatically created and deleted after use. Use this argument
#' to avoid reading files directly from shared locations (namely LabShare). Defaults to FALSE.
#' @return An R list of 16 data frames, one for each chromosome.
#' @examples
#' readall_tab("/Path/to/wiggles/folder/Red1_WT_reps_SacCer3_2mis_MACS_wiggle_norm/")
#' readall_tab("/Path/to/wiggles/", use_readr = TRUE, progress_bar = TRUE)
#' readall_tab("/Path/to/wiggles/", use_readr = T, progress_bar = T, local_copy = T)
#' @export

readall_tab <- function(fileLocation, use_readr = FALSE,
                        progress_bar = FALSE, local_copy = FALSE) {
  ptm <- proc.time()
  if (local_copy) {
    cat('Copying data files to local folder "/temp"\n...')
    # Create temporary directory in current working directory
    # and make it the destination
    dir.create(paste0(getwd(), '/temp'))
    # Copy the files to the new temporary directory
    fileLocation <- '/Volumes/LabShare/HTGenomics/HiSeqOutputs/AveragedReplicates_S288C_SacCer3/Red1_WT_reps_S288C_MACS_wiggle_norm'
    file.copy(fileLocation, 'temp', recursive = TRUE)
    # Update fileLocation to be the local directory
    fileLocation <- paste0('temp/', list.files('temp'))
  }
  
  cat('\nReading data:\n')
  
  filenames <- list.files(fileLocation, full = T)
  if (length(filenames) == 17) {
    filenames <- filenames[2:17]
  }
  if (use_readr) {
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("R package 'readr' required. Please install it:\n", 
           "install.packages('readr')", call. = FALSE)
    }
    library(readr)
    if (progress_bar) {
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("R package 'pbapply' required. Please install it:\n", 
             "install.packages('pbapply')", call. = FALSE)
      }
      library(pbapply)
      alldata <- pblapply(filenames, read_tsv, skip = 2, 
                          col_names = FALSE, progress = FALSE)
    }
    else {
      alldata <- lapply(filenames, read_tsv, skip = 2, 
                        col_names = FALSE, progress = FALSE)
    }
  }
  else {
    if (progress_bar) {
      if (!requireNamespace("pbapply", quietly = TRUE)) {
        stop("R package 'pbapply' required. Please install it:\n", 
             "install.packages('pbapply')", call. = FALSE)
      }
      library(pbapply)
      alldata <- pblapply(filenames, read.table, skip = 2, 
                          sep = "\t")
    }
    else {
      alldata <- lapply(filenames, read.table, skip = 2, 
                        sep = "\t")
    }
  }
  
  if (local_copy) {
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