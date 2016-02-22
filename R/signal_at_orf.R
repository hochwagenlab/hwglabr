#' Signal at all ORFs genome-wide
#'
#' This function allows you to pull out the ChIP signal over all ORFs in the genome. It collects the
#' signal over each ORF plus both flanking regions (1/2 the length of the ORF on each side) and
#' scales them all to the same value: 2 kb. This means that for two example genes with lengths of
#' 500 bp and 2 kb, flanking regions of 250 bp and 1 kb, respectively, will be colleted up and
#' downstream. Both genes will then be scaled to 1 kb and all four flanking regions to 500 bp.
#' The function takes as input the wiggle data as a list of 16 chromosomes (output of 'readall_tab()').
#' @param inputData As a list of the 16 chr wiggle data (output of readall_tab). No default.
#' @param gff Optional dataframe of the gff providing the ORF cordinates. Must be provided if 'gffFile'
#' is not. No default. Note: You can use the function 'gff_read()' in this package to load your selected
#' gff file.
#' @param gffFile Optional string indicating path to the gff file providing the ORF cordinates. Must be
#' provided if 'gff' is not. No default.
#' @param saveFile Boolean indicating whether output should be written to a .txt file (in current working
#' directory). If 'saveFile = FALSE', output is returned to screen or an R object (if assigned).
#' Defaults to FALSE.
#' @return A local data frame with four columns: chr (chromosome number), position (normalized genome
#' between 0 and 2000 bp), signal (ChIP signal) and gene (systematic gene name).
#' @examples
#' signal_at_orf(WT, gff = gff)
#' 
#' signal_at_orf(WT, gffFile = S288C_annotation_modified.gff, saveFile = TRUE)
#' @export

signal_at_orf <- function(inputData, gff, gffFile, saveFile = FALSE) {
  ptm  <- proc.time()
  
  # gff
  if (missing(gffFile) & missing(gff)) {
    stop("No gff data provided.\n",
         "You must provide a lodaded gff as an R data frame ('gff' argument) or the path to a gff file to load ('gffFile' argument).",
         call. = FALSE)
  } else if (!(missing(gffFile) | missing(gff))) {
    stop("Two gff data sources provided.\n",
         "Please provide either a gff R data frame ('gff' argument) or the path to a gff file ('gffFile' argument), not both.\n",
         "Jeez, make up your mind...", call. = FALSE)
  } else if (missing(gff)) {
    gff <- gff_read(gffFile)
    cat('Loaded gff file...\n')
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
    stop("The reference genome in the input data and the gff do not seem to match.\n",
         "Please provide data and gff for the same reference genome.\n", call. = FALSE)
  } else if (check_S288C & check_gff_S288C) {
    cat('Detected ref. genome - S288C\n')
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
  library(dplyr)
  
  cat('Collecting signal...\n')
  cat('(Skip genes whose start/stop coordinates are not found in wiggle data)\n')

  # Create data frames to collect final data for all chrs
  plus_final <- data.frame()
  minus_final <- data.frame()
  
  for(i in 1:length(chrom)) {
    chrNum <- paste0('chr', chrom[i])
    cat(paste0(chrNum, ':\n'))
    
    # Index of ChIP data list item corresponding to chrom to analyze
    # Add '.' to make it unique (otherwise e.g. 'chrI' matches 'chrII' too)
    listIndex <- grep(paste0(chrNum, '.'), names(inputData), fixed = TRUE)
    chromData <- inputData[[listIndex]]
    colnames(chromData) <- c('position', 'signal')
    
    ############################### plus strand ################################
    # Create data frame to collect final data for all genes in chr
    plus_sigs <- data.frame()
    # Get all genes on "+" strand of current chromosome
    chromGff <- gff[gff[, 1] == chrNum & gff[, 7] == '+', ]
    
    geneCount <- 0
    for(j in 1:nrow(chromGff)) {
      # Skip if gene coordinates not in ChIPseq data
      # Comparison below will not work if dplyr was loaded when the wiggle data was loaded
      # (because of dplyr's non standard evaluation: data is in tbl_df class)
      # Workaround: convert to data.frame first
      if(!chromGff[j, 4] %in% as.data.frame(chromData)[, 1] |
         !chromGff[j, 5] %in% as.data.frame(chromData)[, 1]) {
        next
      }
      # Collect flanking regions scaled according to ratio gene length / 1 kb
      gene_leng <- chromGff[j, 5] - chromGff[j, 4]
      start <- chromGff[j, 4] - round((0.5 * gene_leng))
      end <- chromGff[j, 5] + round((0.5 * gene_leng))
      full_leng <- end - start
      gene <- chromGff[j, 9]
      
      # pull out red1 signal
      sig_gene <- chromData[which(chromData[, 1] >= start & chromData[, 1] <= end), ]
      
      # normalize to 2 kb segment length
      sig_gene[, 1] <- sig_gene[, 1] - start
      sig_gene[, 1] <- round(sig_gene[, 1] * (2000 / full_leng))
      
      chr <- as.data.frame(rep(paste0('chr', chrom[i]), nrow(sig_gene)))
      colnames(chr) <- "chr"
      g <- as.data.frame(rep(gene, nrow(sig_gene)))
      colnames(g) <- "gene"
      all <- bind_cols(chr, sig_gene, g)
      # To collect all genes
      plus_sigs <- bind_rows(plus_sigs, all)
      
      geneCount <- geneCount + 1
    }
    cat(paste0('... + strand: ', geneCount, ' genes (skipped ', j - geneCount, ')\n'))
    
    # To collect all chrs
    plus_final <- bind_rows(plus_final, plus_sigs)
    
    ############################## minus strand ##################################
    # Create data frame to collect final data for all genes in chr
    minus_sigs <- data.frame()
    # Get all genes on "-" strand of current chromosome
    chromGff <- gff[gff[, 1] == chrNum & gff[, 7] == '-', ]
    
    geneCount <- 0
    for(j in 1:nrow(chromGff)) {
      # Skip if gene coordinates not in ChIPseq data
      # Comparison below will not work if dplyr was loaded when the wiggle data was loaded
      # (because of dplyr's non standard evaluation: data is in tbl_df class)
      # Workaround: convert to data.frame first
      if(!chromGff[j, 4] %in% as.data.frame(chromData)[, 1] |
         !chromGff[j, 5] %in% as.data.frame(chromData)[, 1]) {
        next
      }    
      # Collect flanking regions scaled according to ratio gene length / 1 kb
      gene_leng = chromGff[j, 5] - chromGff[j, 4]
      start <- chromGff[j, 4] - round((0.5 * gene_leng))
      end <- chromGff[j, 5] + round((0.5 * gene_leng))
      full_leng <- end - start
      gene <- chromGff[j, 9]
      
      # Pull out red1 signal
      sig_gene <- chromData[which(chromData[, 1] >= start & chromData[, 1] <= end), ]
      
      # Normalize to 2 kb segment length
      sig_gene[, 1] <- sig_gene[, 1] - start
      sig_gene[, 1] <- round(sig_gene[, 1] * (2000 / full_leng))
      
      # Reverse the order of the position values
      sig_gene[, 1] <- 2000 - sig_gene[, 1]
      
      chr <- as.data.frame(rep(paste0('chr', chrom[i]), nrow(sig_gene)))
      colnames(chr) <- "chr"
      g <- as.data.frame(rep(gene, nrow(sig_gene)))
      colnames(g) <- "gene"
      all <- bind_cols(chr, sig_gene, g)
      # To collect all genes
      minus_sigs <- bind_rows(minus_sigs, all)
      
      geneCount <- geneCount + 1
    }
    # To collect all chrs
    minus_final <- bind_rows(minus_final, minus_sigs)
    
    cat(paste0('... - strand: ', geneCount, ' genes (skipped ', j - geneCount, ')\n'))
  }
  
  # Merge '+' and '-' strand data
  mergedStrands <- bind_rows(plus_final, minus_final)
  colnames(mergedStrands) <- c("chrom", "position", "signal", "gene")
  # Sort by gene and position
  mergedStrands <- mergedStrands[order(mergedStrands$gene, mergedStrands$position), ]
  
  if(saveFile) {
    cat(paste0('Saving file...\n'))
    if(check_S288C) {
      write.table(mergedStrands, paste0(deparse(substitute(inputData)),
                                        "_S288C_ORFsignal.txt"), sep = "\t", quote = FALSE)
    } else {
      write.table(mergedStrands, paste0(deparse(substitute(inputData)),
                                        "_SK1_ORFsignal.txt"), sep = "\t", quote = FALSE)
    }
  } else {
    return(mergedStrands)
  }
  
  cat(paste0('Completed in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.\n'))
}