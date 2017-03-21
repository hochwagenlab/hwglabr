#' Load wiggle data
#'
#' This function allows you to load tab-separated wiggle data.
#' Adapted from function written by Tovah Markowitz (original function name: \code{readAll.tab}).
#' @param fileLocation A string with the path to the folder containing the wiggle data. No default.
#' @param useReadr Boolean indicating whether to use the much faster \code{read_tsv()} from
#' Hadley Wickham's readr package instead of base R's \code{read.table}. Defaults to \code{TRUE}.
#' @param progressBar Boolean indicating whether to display a progress bar (using R package
#'  pbapply). Defaults to \code{TRUE}.
#' @param localCopy Boolean indicating whether to create local copy of target files before reading
#' data. If \code{TRUE} a local folder is automatically created and deleted after use. Use this argument
#' to avoid reading files directly from shared locations (namely LabShare). Defaults to \code{FALSE}.
#' @param asBedGraph Boolean indicating whether to return the data in bedGraph-like format instead of
#' R list of 16 data frames. Genomic coordinates are kept as 1-based. Defaults to \code{FALSE}.
#' @return Either an R list of 16 data frames, one for each chromosome (\code{asBedGraph=FALSE})
#' or a bedGraph-like data frame (\code{asBedGraph=TRUE}).
#' @examples
#' \dontrun{
#' readall_tab("/Path/to/wiggles/folder/Red1_WT_reps_SacCer3_2mis_MACS_wiggle_norm/")
#' readall_tab("/Path/to/wiggles/", useReadr = TRUE, progressBar = TRUE)
#' readall_tab("/Path/to/wiggles/", useReadr = T, progressBar = T, localCopy = T)
#' #' readall_tab("/Path/to/wiggles/", asBedgraph = TRUE)
#' }
#' @export

readall_tab <- function(fileLocation, useReadr = TRUE,
                        progressBar = TRUE, localCopy = FALSE,
                        asBedGraph = FALSE) {
  ptm <- proc.time()
  
  # Check that path to files is correct
  if (!file.exists(fileLocation)) {
    stop('Cannot seem to find the files.\n',
         'Please check that the provided path to files is correct.', call. = FALSE)
  }
  
  if (localCopy) {
    message('Copying data files to local folder "/readall_tab_temp"\n...')
    # Check if a directory 'readall_tab_temp' already exists
    if (file.exists('readall_tab_temp')) {
      stop('A folder called "readall_tab_temp" already exists in the current working directory.\n',
           'Please remove it and repeat function call.', call. = FALSE)
    }
    # Create temporary directory in current working directory
    # and make it the destination
    dir.create('readall_tab_temp')
    # Copy the files to the new temporary directory
    file.copy(fileLocation, 'readall_tab_temp', recursive = TRUE)
    # Update fileLocation to be the local directory
    fileLocation <- paste0('readall_tab_temp/', list.files('readall_tab_temp'))
  }
  
  message('Reading data...')
  
  filenames <- list.files(fileLocation, full.names = T)
  if (length(filenames) == 17) {
    filenames <- filenames[2:17]
  } else stop('The provided "fileLocation" argument seems to be wrong: found ', length(filenames), 
         ' files. Please provide the path to the folder containing the wiggle files.', call. = FALSE)
  
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
      # Supress messages in order to suppress readr's printed col specs
      suppressMessages(alldata <- pbapply::pblapply(filenames, readr::read_tsv, skip = 2,
                                                    col_names = FALSE, progress = FALSE))
    }
    else {
      # Supress messages in order to suppress readr's printed col specs
      suppressMessages(alldata <- lapply(filenames, readr::read_tsv, skip = 2,
                                         col_names = FALSE, progress = FALSE))
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
    message('Deleting local copy\n...')
    # Delete temporary local directory
    unlink('readall_tab_temp', recursive = TRUE)
  }
  
  # Check which reference genome was used to map seq. data
  names(alldata) = filenames
  check_S288C <- any(grep("chrI.", names(alldata), fixed = TRUE))
  check_SK1 <- any(grep("chr01.", names(alldata), fixed = TRUE))
  if (check_S288C) {
    message("Detected ref. genome - S288C\n(Chrs numbered using roman numerals)")
  }
  else if (check_SK1) {
    message("Detected ref. genome - SK1\n(Chrs numbered using arabic numerals)")
  }
  else stop("Did not recognize reference genome.")
  
  # Convert to bedGraph-like format
  message('Converting to bedGraph-like...')
  if(asBedGraph){
    get_chr <- function(list_element_name, regex="chr[IXV]*"){
      index <- gregexpr(regex, list_element_name)[[1]]
      length <- attributes(index)$match.length - 1
      return(substr(list_element_name, index[1], index + length))
    }
    
    if (check_S288C) {
      regex <- "chr[IXV]*"
    }
    else if (check_SK1) {
      regex <- "chr[[:digit:]]*"
    }
    
    for(i in 1:length(alldata)){
      alldata[[i]] <- data.frame(get_chr(names(alldata)[i], regex=regex),
                                 alldata[[i]][, 1], alldata[[i]][, 1] + 1,
                                 alldata[[i]][, 2], row.names = NULL)
    }
    
    # Collapse to df (use data.table package if available)
    if(!requireNamespace("data.table", quietly = TRUE)){
      message("...")
      message("Note:")
      message("Install package 'data.table' to significantly decrease this function's runtime.")
      message("...")
      alldata <- do.call('rbind', alldata)
    } else alldata <- data.table::rbindlist(alldata)
  }
  
  colnames(alldata) <- c('chr', 'start', 'end', 'score')
  
  elapsed_time <- proc.time()[3] - ptm[3]
  
  # Convert to appropriate unit
  if (elapsed_time < 60) {
    elapsed_time <- paste0(round(elapsed_time, 1), " sec.")
  } else if (elapsed_time >= 60 & elapsed_time < 3600) {
    elapsed_time <- paste0(round(elapsed_time / 60, 1), " min.")
  } else {
    elapsed_time <- paste0(round(elapsed_time / 60 / 60, 1), " h.")
  }
  
  message("\n...\nCompleted in ", elapsed_time)
  return(alldata)
}
