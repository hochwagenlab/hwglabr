#' Standard analysis of HT genomics experiments
#'
#' This function will run the lab's standard analysis of HT genomics experiments for which
#' tab-separated wiggle data is generated. It will call different functions in the package
#' to produce several .pdf files of analysis plots, written to a new folder in
#' ".../LabShare/HTGenomics/Opening_act/".
#' @param wiggleData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}).
#' No default.
#' @param relevantGenotype String indicating the relevant strain mutations. Just use "WT"
#' if there are no relevant mutations. No default.
#' @param chipTarget String indicating the ChIP target protein. No default.
#' No default.
#' @param sampleID String indicating the sample ID, including the ID used in the
#' analysis pipeline (with a date) and the read mapping conditions (see examples below).
#' The function asks the user to check that the provided "sampleID" matches the required
#' format before proceeding with the analysis. No default.
#' @param userInput Boolean indicating whether to ask user to check the format of the
#' \code{sampleID} argument. Defaults to \code{TRUE}.
#' No default.
#' @return A new folder in ".../LabShare/HTGenomics/Opening_act/" containing several
#' plots as .pdf files.
#' @examples
#' \dontrun{
#' opening_act(wiggleData=WT, relevantGenotype="WT", chipTarget="Red1",
#'             sampleID="AH119C-040114-sacCer3-2mis")
#' opening_act(set1_wiggle_data, "set1", "Red1", "AH8584b-16032016-sacCer3-2mis")
#' opening_act(rec8, "rec8", "Red1", "AH8115b-24042015-SacCer3-2mis", userInput=FALSE)
#' }
#' @export

opening_act <- function(wiggleData, relevantGenotype, chipTarget, sampleID,
                        userInput = TRUE) {
  ptm <- proc.time()
  
  if(userInput){
    # Ask user to make sure they provided a valid ID for the data set
    title <- paste0('The "sampleID" argument will be used to name the final output folder.
It should identify the yeast strain, date, and read mapping conditions, as in:
"AH119C-040114-sacCer3-2mis".\n
You provided the string "', sampleID, '" as the sampleID. Is this correct?')
    choices = c('No, let me change that.',
                'Yes, continue analysis!')
    answer <- menu(choices, graphics = FALSE, title)
    
    if(answer == 0 | answer == 1){
      stop('You chose to stop the function.', call. = FALSE)
    }
  }
  
  # Check which reference genome was used to map seq. data
  check_S288C <- any(grep('chrI.', names(wiggleData), fixed = TRUE))
  check_SK1 <- any(grep('chr01.', names(wiggleData), fixed = TRUE))
  if (check_S288C) {
    refGenome <- 'S288c'
    message("Detected ref. genome - S288c (Chrs numbered using roman numerals)")
    
    if(exists('s288C_gff')){
      gff_file <- s288C_gff
    } else {
      stop('Cannot load data. Please load the package before running:
           library(hwglabr)', call. = FALSE)
    }
  }
  else if (check_SK1) {
    refGenome <- 'SK1'
    message("Detected ref. genome - SK1 (Chrs numbered using arabic numerals)")
    
    if(exists('SK1_gff')){
      gff_file <- SK1_gff
    } else {
      stop('Cannot load data. Please load the package before running:
           library(hwglabr)', call. = FALSE)
    }
  }
  else stop("Did not recognize reference genome.
            Please make sure chromosome numbers follow the standard format.",
            call. = FALSE)
  
  # Check that gff data (for signal_at_orf) can be loaded
  
  
  
  destination <- "/Volumes/LabShare/HTGenomics/Opening_act/"
  output_dir <- paste0(relevantGenotype, '_anti-', chipTarget, '_',
                       sampleID)
  # Check if the directory already exists
  if (file.exists(paste0(destination, output_dir))) {
    stop('A folder named "', output_dir, '" already exists in "Opening_act".\n',
         call. = FALSE)
  }
  # Create output directory
  message('Creating output directory "', output_dir, '"')
  dir.create(file.path(paste0(destination, output_dir)))
  
  #----------------------------------------------------------------------------#
  #                                Run analysis                                #
  #----------------------------------------------------------------------------#
  
  #----------------------------------------------------------------------------#
  # Chr size bias
  message('... Chromosome size bias:')
  suppressMessages(output <- hwglabr::chr_coverage(wiggleData))

  suppressMessages(
    hwglabr::chr_coverage_plot(output, genome = refGenome, onScreen = FALSE,
                               fileName = paste0(destination, output_dir, '/',
                                                 output_dir, '_sizeBias.pdf'))
  )
  
  message('Saved plot ', paste0(output_dir, '_sizeBias.pdf'))
  
  
  #----------------------------------------------------------------------------#
  # Signal at centromere
  message('... Signal at centromeres:')
  
  # convert S288Ccen or SK1cen into bed file
  if (check_S288C) {
    cen <- S288Ccen
    
  } else if (check_SK1) {
    cen <- SK1cen

  } else stop("Did not recognize reference genome.")
  
  cenBed <- data.frame(cen$Chromosome, cen$Mid, cen$Mid + 1, stringsAsFactors=F)
  
  # calculate average around centromere
  suppressMessages(
    wiggle_cen_avg <- signal_average( signal_at_summit(wiggleData, cenBed, 50000,
                                                       onlyComplete=F) )
  )
  
  # plot results
  fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtCen.pdf')
  pdf(file = paste0(fileName), width = 8, height = 4)
  
  YLIM <- range(wiggle_cen_avg$mean_signal)
  if( YLIM[[2]] < 2) { YLIM[2] <- 2 }
  plot(wiggle_cen_avg$position/1000, wiggle_cen_avg$mean_signal, type="l",
       ylim=YLIM, xlab="Distance around Centromere (kb)", ylab="Signal", 
       lwd=1, cex.axis=1, las=1, col="darkorange", cex.lab=1,
       main=paste0("Signal around centromeres: ", refGenome), cex.main=1)
  
  dev.off()
  message('Saved plot ', paste0(output_dir, '_signalAtCen.pdf'))
  
  #----------------------------------------------------------------------------#
  # Signal at rDNA
  message('... Signal flanking rDNA:')
  suppressMessages(rDNA <- signal_at_rDNA(wiggleData, saveFile = F))
  colnames(rDNA) <- c('position', 'signal')
  
  # plot results
  fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtrDNA.pdf')
  pdf(file = paste0(fileName), width = 6, height = 3)
  
  if (check_S288C) {
    start <- 451575
    end <- 468931
  } else {
    start <- 433029
    end <- 451212
  }
  
  plot(rDNA$position/1000, rDNA$signal, type="l",
       xlab="Position on chr 12 (kb)", ylab="Signal", 
       lwd=1, cex.axis=1, las=1, col='black', cex.lab=1,
       main=paste0("Signal around rDNA: ", refGenome), cex.main=1)
  # Add labels for rDNA
  axis(1, side=3, at = c(start / 1000, end / 1000),
       labels = c('rDNA >', '< rDNA'),
       col = 'red', lwd = 3)
  
  dev.off()
  message('Saved plot ', paste0(output_dir, '_signalAtrDNA.pdf'))
  
  #----------------------------------------------------------------------------#
  # Signal from telomeres
  message('... Signal at sub-telomeric regions:')
  
  
  
  #----------------------------------------------------------------------------#
  # Jonna's
  message('... Signal at :')
  
  
  #----------------------------------------------------------------------------#
  # Meta ORF
  message('... Signal at meta ORF analysis:')
  meta_orf <- hwglabr::signal_at_orf(wiggleData, gff = gff_file, saveFile = F)
  suppressMessages(meta_orf <- hwglabr::signal_average(meta_orf, saveFile = F))
  
  # plot results
  fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtORF.pdf')
  pdf(file = paste0(fileName), width = 4, height = 3)
       
  plot(meta_orf, type = 'l', xaxt = 'n', yaxt = 'n',
       xlim = c(0, 1000), lwd = 4, col = 'orange',
       xlab = "Scaled ORF",
       ylab = 'Signal', main = paste0('Signal at meta ORF'), bty = "n")
  axis(1, at = c(0, 250, 750, 1000), labels = c('', 'start', 'stop', ''), las = 1)
  axis(2, las = 2)
  abline(v = c(250, 750), lty= 2)
  dev.off()
  message('Saved plot ', paste0(output_dir, '_signalAtORF.pdf'))
  
  #----------------------------------------------------------------------------#
  
  message()
  message('------------------')
  message('All plots saved to ', paste0(destination, output_dir))
  message('------------------')
  
  elapsed_time <- round((proc.time()[3] - ptm[3]), 1)
  if(elapsed_time < 60){
    message('\n...\nCompleted in ', elapsed_time, ' sec.')
  } else message('\n...\nCompleted in ', elapsed_time, ' min.')
}