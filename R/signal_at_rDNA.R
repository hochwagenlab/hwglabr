#' Signal flanking rDNA
#'
#' This function allows you to pull out the ChIP signal flanking the rDNA region on
#' chromosome 12. The function takes as input the wiggle data as a list of 16 chromosomes.
#' (output of \code{\link{readall_tab}}).\cr
#' @param inputData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current working
#' directory). If \code{saveFile = FALSE}, output is returned to screen or an R object (if assigned).
#' Defaults to \code{FALSE}.
#' @return A local data frame with four columns:
#' \enumerate{
#'   \item \code{chr} Chromosome number
#'   \item \code{position} Nucleotide coordinate (in normalized total length of 1 kb)
#'   \item \code{signal} ChIP-seq signal at each position (1 to 1000)
#'   \item \code{gene} Systematic gene name
#' }
#' @examples
#' signal_at_rDNA(WT)
#' 
#' signal_at_rDNA(WT, saveFile = TRUE)
#' @export

signal_at_rDNA <- function(inputData, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # Check reference genome
  chrom_S288C <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                   "XI", "XII", "XIII", "XIV", "XV", "XVI")
  chrom_SK1 <- c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                 '11', '12', '13', '14', '15', '16')
  
  check_S288C <- any(grep('chrI.', names(inputData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(inputData), fixed = TRUE))
  
  if (check_S288C) {
    refGenome <- 'S288c'
    message("Detected ref. genome - S288c (Chrs numbered using roman numerals)")
  } else if (check_SK1) {
    refGenome <- 'SK1'
    message("Detected ref. genome - SK1 (Chrs numbered using arabic numerals)")
  } else stop("Did not recognize reference genome.
              Please make sure chromosome numbers follow the standard format, e.g. 'chrI' or 'chr01'.",
              call. = FALSE)
  
  message('\nCollecting signal...')
  
  ### rDNA coordinates collected from gff files:
  # Collect signal +/- 40 kb flanking rDNA and get index of list element containing chr12
  if (check_S288C) {
    start <- 451575 - 40000
    end <- 468931 + 40000
    chr12 <- grep('chrXII.', names(inputData), fixed = TRUE)
  } else {
    start <- 433029 - 40000
    end <- 451212 + 40000
    chr12 <- grep('chr01.', names(inputData), fixed = TRUE)
  }
  
  # Keep only chr12
  inputData <- inputData[[chr12]]
  
  # Convert positions to indices
  start <- tail(inputData[inputData[, 1] < start, 1], 1)
  start <- which(inputData[, 1] == start[[1]])
  
  end <- head(inputData[inputData[, 1] > end, 1], 1)
  end <- which(inputData[, 1] == end[[1]])
  
  # Get subset of data corresponding to rDNA +/- 40 kb
  rDNA_signal <- inputData[start:end, ]
  colnames(rDNA_signal) <- c('position_on_chr12', 'signal')

  
  if(saveFile) {
    message(paste0('Saving file...\n'))
    if(check_S288C) {
      write.table(rDNA_signal,
                  paste0(deparse(substitute(inputData)), "_S288C_signalAtrDNA.txt"),
                  sep = "\t", quote = FALSE,
                  row.names = FALSE)
      message('Done!')
    } else {
      write.table(rDNA_signal,
                  paste0(deparse(substitute(inputData)), "_S288C_signalAtrDNA.txt"),
                  sep = "\t", quote = FALSE,
                  row.names = FALSE)
      message('Done!')
    }
  } else {
    message('Done!')
    return(rDNA_signal)
  }
}
