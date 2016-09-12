#' Mean signal for all ORFs genome-wide
#'
#' This function allows you to calculate mean ChIP signal for each ORF in the genome. The function
#' goes through all included features in the supplied gff data and for each one it adds ORF length
#' (in bp) and the mean of the signal collected from the supplied ChIP-seq data. This can then be
#' used to analyse mean signal as a function of ORF length.\cr
#' The function takes as input the wiggle data as a list of 16 chromosomes (output of
#' \code{\link{readall_tab}}).
#' \cr \cr
#' \strong{Note:} Our wiggle data always contains gaps with missing chromosome coordinates
#' and ChIP-seq signal. The way this function deals with that is by skipping affected genes.
#' The number of skipped genes in each chromosome is printed to the console, as well as the
#' final count (and percentage) of skipped genes. \cr
#' @param inputData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.
#' @param gff Optional dataframe of the gff providing the ORF cordinates. Must be provided if
#' \code{gffFile} is not. No default. Note: You can use the function \code{\link{gff_read}} in hwglabr to
#' load your selected gff file.
#' @param gffFile Optional string indicating path to the gff file providing the ORF cordinates. Must be
#' provided if \code{gff} is not. No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current working
#' directory). If \code{saveFile = FALSE}, output is returned to screen or an R object (if assigned).
#' Defaults to FALSE.
#' @return A data frame equivalent to the supplied gff table with the following additional columns:
#' \enumerate{
#'   \item \code{length} Length of the ORF in bp
#'   \item \code{mean_signal} Mean of the signal collected from the supplied data
#' }
#' \strong{Note:} Skipped genes are included in the output with 'NA' for \code{mean_signal}.
#' @examples
#' \dontrun{
#' signal_per_orf_length(WT, gff = gff)
#' 
#' signal_per_orf_length(WT, gffFile = S288C_annotation_modified.gff, saveFile = TRUE)
#' }
#' @export

signal_per_orf_length <- function(inputData, gff, gffFile, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # gff
  if (missing(gffFile) & missing(gff)) {
    stop("No gff data provided.\n",
         "You must provide either a lodaded gff as an R data frame ('gff' argument)
         or the path to a gff file to load ('gffFile' argument).",
         call. = FALSE)
  } else if (!(missing(gffFile) | missing(gff))) {
    stop("Two gff data sources provided.\n",
         "Please provide either a gff R data frame ('gff' argument)
         or the path to a gff file ('gffFile' argument), not both.\n",
         call. = FALSE)
  } else if (missing(gff)) {
    gff <- hwglabr::gff_read(gffFile)
    message('Loaded gff file...\n')
  }
  
  # Check reference genome for both the input data and the gff file; make sure they match
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  
  check_S288C <- any(grep('chrI.', names(inputData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(inputData), fixed = TRUE))
  
  check_gff_S288C <- any(gff[, 1] == 'chrI')
  check_gff_SK1 <- any(gff[, 1] == 'chr01')
  
  if (check_S288C != check_gff_S288C) {
    stop("The reference genomes in the input data and the gff do not seem to match.\n",
         "Please provide data and gff for the same reference genome.\n", call. = FALSE)
  } else if (check_S288C & check_gff_S288C) {
    message('Detected ref. genome - S288C\n')
    chrom <- chrom_S288C
  } else if (check_SK1 & check_gff_SK1) {
    print('Detected ref. genome - SK1')
    chrom <- chrom_SK1
  } else stop("Did not recognize reference genome.
              Check that chromosome numbers are in the usual format, e.g. 'chrI' or 'chr01'.")
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' needed for this function to work. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  message('\nThe following types of features are present in the gff data you provided
      (they will all be included in the analysis):\n')
  for(i in 1:length(unique(gff[, 3]))) {
    message(unique(gff[, 3])[i])
  }
  
  message('\nCollecting signal...')
  message('Skip ORFs with missing coordinates and signal in wiggle data)\n')
  
  # Add gene length to gff plus a new column for the mean signal
  gff$length <- gff$end - gff$start
  gff$mean_signal <- NA
  
  # Create data frame to collect final data for all chrs
  gff_final <- data.frame()
  # Keep track of total and non-skipped genes, to print info at the end
  number_genes <- 0
  number_skipped_genes <- 0
  
  # Iterate over chrs
  for(i in 1:length(inputData)) {
    chrNum <- paste0('chr', chrom[i])
    message(paste0(chrNum, ':\n'))
    
    # Index of ChIP data list item corresponding to chrom to analyze
    # Add '.' to make it unique (otherwise e.g. 'chrI' matches 'chrII' too)
    listIndex <- grep(paste0(chrNum, '.'), names(inputData), fixed = TRUE)
    chromData <- inputData[[listIndex]]
    colnames(chromData) <- c('position', 'signal')
    
    # Get all genes on current chromosome
    chromGff <- gff[gff[, 1] == chrNum, ]
    
    # Count skipped genes
    geneCount <- 0
    for (j in 1:nrow(chromGff)) {
      # Skip if gene coordinates not in ChIPseq data
      # Comparison below will not work if dplyr was loaded when the wiggle data was loaded
      # (because of dplyr's non standard evaluation: data is in tbl_df class)
      # Workaround: convert to data.frame first
      if(!chromGff[j, 4] %in% as.data.frame(chromData)[, 1] |
         !chromGff[j, 5] %in% as.data.frame(chromData)[, 1]) {
        next
      }
      
      # Collect signal
      signal_vector <- as.data.frame(chromData[chromData[, 1] >= chromGff$start[j] &
                                   chromData[, 1] <= chromGff$end[j], 2])
      
      # Skip if there are discontinuities in the data (missing position:value pairs)
      if(nrow(signal_vector) != chromGff$length[j] + 1) next
      
      # Calculate mean signal on gene and save in table
      chromGff$mean_signal[j] <- mean(signal_vector[, 1])
      
      geneCount <- geneCount + 1
    }
    
    # To collect all chrs
    gff_final <- dplyr::bind_rows(gff_final, chromGff)
    
    # Keep track of total and non-skipped genes, to print info at the end
    number_genes <- number_genes + j
    number_skipped_genes <- number_skipped_genes + (j - geneCount)
    
    message(paste0('... ', geneCount, ' ORFs (skipped ', j - geneCount, ')\n'))
  }
  # Print info on total and non-skipped genes
  message('------')
  message(paste0('Skipped ', number_skipped_genes, ' of a total of ', number_genes,
             " ORFs (", round((number_skipped_genes * 100 / number_genes), 1),
             "%).\n"))
  message('------')
  
  message(paste0('Completed in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.\n'))
  
  if(saveFile) {
    message(paste0('Saving file...\n'))
    if(check_S288C) {
      write.table(gff_final, paste0(deparse(substitute(inputData)),
                                        "_S288C_mean_signal_perORF.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      message('Done!')
    } else {
      write.table(mergedStrands, paste0(deparse(substitute(inputData)),
                                        "_SK1__mean_signal_perORF.txt"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      message('Done!')
    }
  } else {
    message('Done!')
    return(gff_final)
  }
}