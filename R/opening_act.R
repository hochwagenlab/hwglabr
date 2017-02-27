#' Standard analysis of ChIP-seq experiments
#'
#' This function will run the lab's standard analysis of ChIP-seq experiments, for which
#' tab-separated wiggle data is generated. It will call different functions in the package
#' to produce several .pdf files of analysis plots, written to a new folder in
#' ".../LabShare/HTGenomics/Opening_act/".
#' \cr \cr
#' \strong{Note:} When running on data aligned to the SK1 genome three of the plots are
#' not produced: signal at sub-telomeric regions, signal at DSB hotspots and signal at
#' axis binding sites. The reason for skipping the sub-telomeric signal analysis is the
#' fact that the SK1 genome annotation contains inconsistencies at sub-telomeric regions.
#' The reason for skipping the other two analysis is the fact that the reference data for
#' DSB hotspots and Red1 binding sites were not available for SK1 at the time of writing
#' this function. It is probably redundant to run this (long!) analysis for both genomes
#' anyway, and using data mapped to the S288c reference genome should be preferred.\cr
#' @param wiggleData As a list of the 16 chr wiggle data (output of \code{\link{readall_tab}}).
#' No default.
#' @param relevantGenotype String indicating the relevant strain mutations. Just use "WT",
#' for example, if there are no relevant mutations. No default.
#' @param chipTarget String indicating the ChIP target protein. No default.
#' @param sampleID String indicating the sample ID, including the ID used in the
#' analysis pipeline (with a date) and the read mapping conditions (see examples below).
#' The function asks the user to check that the provided "sampleID" matches the required
#' format before proceeding with the analysis. No default.
#' @param userInput Boolean indicating whether to ask user to check the format of the
#' \code{sampleID} argument. Defaults to \code{TRUE}.
#' @param runMetaORF Boolean indicating whether to run the meta ORF analysis. This analysis
#' typically takes about 30 minutes to run, so it may be useful to exclude it.
#' Defaults to \code{TRUE}.
#' @return A new folder in ".../LabShare/HTGenomics/Opening_act/" containing output plots
#' (as .pdf files) of the following analysis:
#' \enumerate{
#'   \item \strong{Chromosome size bias}
#'   \item \strong{Signal at centromeres}
#'   \item \strong{Signal flanking rDNA}
#'   \item \strong{Signal at sub-telomeric regions} (data mapped to S288c reference genome only)
#'   \item \strong{Signal at DSB hotspots} (data mapped to S288c reference genome only)
#'   \item \strong{Signal at axis binding sites} (data mapped to S288c reference genome only)
#'   \item \strong{Signal at meta ORF}
#' }
#' @examples
#' \dontrun{
#' opening_act(wiggleData=WT, relevantGenotype="WT", chipTarget="Red1",
#'             sampleID="AH119C-040114-sacCer3-2mis")
#' opening_act(set1_wiggle_data, "set1", "Red1", "AH8584b-16032016-sacCer3-2mis")
#' opening_act(rec8, "rec8", "Red1", "AH8115b-24042015-SacCer3-2mis",
#'             userInput=FALSE, runMetaORF=FALSE)
#' }
#' @export

opening_act <- function(wiggleData, relevantGenotype, chipTarget, sampleID,
                        userInput = TRUE, runMetaORF = TRUE) {
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
  suppressMessages(output <- hwglabr::chr_coverage(wiggleData, meanNorm = T))

  suppressMessages(
    hwglabr::chr_coverage_plot(output, genome = refGenome, meanNorm = TRUE,
                               onScreen = FALSE,
                               fileName = paste0(destination, output_dir, '/',
                                                 output_dir, '_chrSizeBias.pdf'))
  )
  
  message('Saved plot ', paste0(output_dir, '_chrSizeBias.pdf'))
  
  
  #----------------------------------------------------------------------------#
  # Signal at centromere (Tovah)
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
  pdf(file = paste0(fileName), width = 6, height = 3)
  
  YLIM <- range(wiggle_cen_avg$mean_signal)
  if( YLIM[[2]] < 2) { YLIM[2] <- 2 }
  plot(wiggle_cen_avg$position/1000, wiggle_cen_avg$mean_signal, type="l",
       ylim=YLIM, xlab="Distance around Centromere (kb)", ylab="Signal", 
       lwd=1, cex.axis=1, las=1, col="darkorange", cex.lab=1,
       main=paste0("Signal around centromeres: ", refGenome), cex.main=1)
  
  dev.off()
  message('    Saved plot ', paste0(output_dir, '_signalAtCen.pdf'))
  
  #----------------------------------------------------------------------------#
  # Signal at rDNA
  message('... Signal flanking rDNA:')
  suppressMessages(rDNA <- hwglabr::signal_at_rDNA(wiggleData, saveFile = F))
  colnames(rDNA) <- c('position', 'signal')
  
  # plot results
  fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtrDNA.pdf')
  pdf(file = paste0(fileName), width = 6, height = 3)
  
  plot(rDNA$position/1000, rDNA$signal, type="l",
       xlab="Position on chr 12 (kb)", ylab="Signal", 
       lwd=1, cex.axis=1, las=1, col='black', cex.lab=1, cex.main=1)
  
  # A stretch present in S288C downstream of rDNA is absent in SK1 (the strain we use)
  # Add label for that (from end of rDNA until about bp 490'500)
  if (check_S288C) {
    start <- 468931
    end <- 490500
    axis(1, at = c(start / 1000, end / 1000),
         labels = c('', ''),
         col = 'blue', lwd = 3)
  }
  
  # Add labels for rDNA
  if (check_S288C) {
    start <- 451575
    end <- 468931
    title(paste0("Signal around rDNA: ", refGenome, '\n(rDNA position marked in red;',
                 '\nregion absent form SK1 genome marked in blue)'))
  } else {
    start <- 433029
    end <- 451212
    title(paste0("Signal around rDNA: ", refGenome, '\n(rDNA position marked in red)'))
  }
  
  axis(1, at = c(start / 1000, end / 1000),
       labels = c('', ''),
       col = 'red', lwd = 3)
      
  
  dev.off()
  message('    Saved plot ', paste0(output_dir, '_signalAtrDNA.pdf'))
  
  #----------------------------------------------------------------------------#
  # Signal from telomeres (Viji)
  if (check_S288C) {
    message('... Signal at sub-telomeric regions:')
    
    # Call signal_from_telomeres() function
    suppressMessages(sample_telo <- hwglabr::signal_from_telomeres(wiggleData,
                                                                   lengthToCollect=120000))
    
    # Combine data from large and small chromosomes
    data <- dplyr::summarise(dplyr::group_by(do.call('rbind', c(sample_telo$small_chrs,
                                                                sample_telo$large_chrs)),
                                             distance_to_telomere),
                             mean_signal=mean(signal, na.rm = TRUE))
    
    # Calculate genome average of the ChIP signal
    sums <- vector(length=16)
    counts <- vector(length=16)
    for (i in c(1:16)) {
      sums[i] <- sum(wiggleData[[i]][,2])
      counts[i] <- nrow(wiggleData[[i]])
    }
    wiggleDataGenomeAvg <- sum(sums)/sum(counts)
    
    # Subtract genome average signal from each datapoint to normalize to genome average
    data$mean_signal <- data$mean_signal - wiggleDataGenomeAvg
    
    # Smooth data over 25kb regions?
    data <- ksmooth(data$distance_to_telomere, data$mean_signal, bandwidth = 25000)
    averageSubtelomericSignal <- data.frame('distance_from_telomere' = data[[1]],
                                            'signal' = data[[2]])
    
    # Plot results
    fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtTelomeres.pdf')
    pdf(file = paste0(fileName), width = 6, height = 4)
    
    plot(averageSubtelomericSignal$distance_from_telomere / 1000,
         averageSubtelomericSignal$signal, type="l", lwd=2, col='plum4',
         xlab="Distance from telomeres (Kb)", ylab = "Average Enrichment",
         main=paste0("Signal at sub-telomeric regions: ", refGenome), cex.main=1)
    abline(h = 0, lty=3, lwd=1.5)
    dev.off()
    message('    Saved plot ', paste0(output_dir, '_signalAtTelomeres.pdf'))
  } else {
    message('... Skip signal at sub-telomeric regions')
    message('    (Low quality SK1 genome annotation at telomeric and sub-telomeric regions)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal around DSDs by DSB hotspot hotness (Jonna)
  if (check_S288C) {
    message('... Signal at DSB hotspots:')
    
    Spo11_DSBs$V1 <- as.character(Spo11_DSBs$V1)
    # Order by signal to make 8 groups based on hotspot hotness
    sporder <- Spo11_DSBs[order(Spo11_DSBs$V4),]
    
    spo1 <- sporder[1:450,]
    spo2 <- sporder[451:900,]
    spo3 <- sporder[901:1350,]
    spo4 <- sporder[1351:1800,]
    spo5 <- sporder[1801:2250,]
    spo6 <- sporder[2251:2700,]
    spo7 <- sporder[2701:3150,]
    spo8 <- sporder[3151:3599,]
    
    data <- list(spo1, spo2, spo3, spo4, spo5, spo6, spo7, spo8)
    
    message('    Computing signal around DSB hotspots; this takes a few minutes...')
    suppressMessages(data <- lapply(data, function(x) signal_at_summit(wiggleData, x, 1000)))
    suppressMessages(data <- lapply(data, signal_average))
    
    min_data <- sapply(data, function(x) min(x[, 2]))
    max_data <- sapply(data, function(x) max(x[, 2]))
    
    colors <- c("lightblue1", "cadetblue1", "deepskyblue", "deepskyblue3",
                "royalblue", "blue", "blue4", "black")
    
    fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtDSBhotspots.pdf')
    pdf(file = paste0(fileName), width = 6, height = 5)
    
    plot(0, type="l", lwd=3, xlim = c(-1000, 1000), ylab="ChIP-seq signal",
         ylim=c(min(min_data), max(max_data)),
         xlab="Distance from Hotspot Midpoints (bp)")
    for(i in 1:length(data)){
      lines(data[[i]], lwd=3, col=colors[i])
    }
    
    legend("topright", lty=c(1,1), lwd=3, title="Hotspot strength",
           legend=c("weakest", "", "", "", "", "", "", "hottest"),
           bg = "white", col=colors)
    dev.off()
    
    message('    Saved plot ', paste0(output_dir, '_signalAtDSBhotspots.pdf'))
    
  } else {
    message('... Skip signal at DSB hotspots')
    message('    (hotspots only available for data mapped to S288C genome)')
  }
  
  #----------------------------------------------------------------------------#
  # Signal at axis binding sites (Jonna)
  if (check_S288C) {
    message('... Signal at axis binding sites:')
    
    Red1_summits <- Red1_summits[, c(1, 2, 3, 5)]
    Red1_summits[, 1] <- as.character(Red1_summits[,1])
    # Order by signal to make groups based on binding level
    axisorder <- Red1_summits[order(Red1_summits[,4]), ]
    n <- nrow(Red1_summits)/4
    
    axis1 <- axisorder[1:round(n),]
    axis2 <- axisorder[(round(n)+1):round(2*n),]
    axis3 <- axisorder[(round(2*n)+1):round(3*n),]
    axis4 <- axisorder[(round(3*n)+1):round(4*n),]
    
    data <- list(axis1, axis2, axis3, axis4)
    message('    Computing signal around axis binding sites; this takes a few minutes...')
    suppressMessages(data <- lapply(data, function(x) signal_at_summit(wiggleData, x, 1000)))
    suppressMessages(data <- lapply(data, signal_average))
    
    min_data <- sapply(data, function(x) min(x[, 2]))
    max_data <- sapply(data, function(x) max(x[, 2]))
    
    colors <- c("darkolivegreen3", "green3", "darkgreen", "black")
    
    fileName <- paste0(destination, output_dir, '/', output_dir, '_signalAtAxisSites.pdf')
    pdf(file = paste0(fileName), width = 6, height = 5)
    
    plot(0, type="l", lwd=3, ylab="ChIP-seq signal", xlim = c(-1000, 1000),
         ylim=c(min(min_data), max(max_data)),
         xlab="Distance from Axis Midpoints (bp)",
         main = "Signal around axis sites\n(Red1 peaks in WT)")
    
    for(i in 1:length(data)){
      lines(data[[i]], lwd=3, col=colors[i])
    }
    
    legend("topright", lty=c(1, 1), lwd=3, title="Axis binding",
           legend=c("least", "", "", "most"), bg='white', col=colors)
    dev.off()
    
    message('    Saved plot ', paste0(output_dir, '_signalAtAxisSites.pdf'))
    
  } else {
    message('... Skip signal at axis binding sites')
    message('    (binding sites only available for data mapped to S288C genome)')
  }
  
  #----------------------------------------------------------------------------#
  # Meta ORF
  if(runMetaORF){
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
    message('    Saved plot ', paste0(output_dir, '_signalAtORF.pdf')) 
  }
  
  #----------------------------------------------------------------------------#
  
  message()
  message('------------------')
  message('All plots saved to ', paste0(destination, output_dir))
  message('------------------')
  
  elapsed_time <- round((proc.time()[3] - ptm[3]), 1)
  if(elapsed_time < 60){
    message('\n...\nCompleted in ', elapsed_time, ' sec.')
  } else if(elapsed_time >= 60 & elapsed_time < 3600){
    message('\n...\nCompleted in ', round(elapsed_time / 60, 1), ' min.') 
  } else message('\n...\nCompleted in ', round(elapsed_time / 60 / 60, 1), ' h.')
}