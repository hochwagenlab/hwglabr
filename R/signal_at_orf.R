#' Signal at all ORFs genome-wide (meta ORF)
#'
#' This function allows you to pull out the ChIP signal over all ORFs in the genome. It collects the
#' signal over each ORF plus both flanking regions (1/2 the length of the ORF on each side) and
#' scales them all to the same value (1000). This means that for two example genes with lengths of
#' 500 bp and 2 kb, flanking regions of 250 bp and 1 kb, respectively, will be collected up and
#' downstream. Both gene lengths will then be scaled to 500 and all four flanking regions to 250.
#' After scaling, a loess model of the signal is built and used to output predictions of the signal
#' at each position between 1 and 1000. \cr
#' The function takes as input the wiggle data as a list of 16 chromosomes.
#' (output of \code{\link{readall_tab}}).
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
#' @param loessSpan Number specifying \code{span} argument for \code{loess} function (the smoothing
#' parameter alpha). This controls the degree of smoothing of the signal.
#' Defaults to \code{0.05}.
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
#' \dontrun{
#' signal_at_orf(WT, gff = gff)
#' 
#' signal_at_orf(WT, gffFile = S288C_annotation_modified.gff,
#'               loessSpan = 0.1, saveFile = TRUE)
#' }
#' @export

signal_at_orf <- function(inputData, gff, gffFile, loessSpan = 0.05, saveFile = FALSE) {
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
    message('Loaded gff file...')
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
    message('Detected ref. genome - S288C')
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
(they will all be included in the analysis):')
  for(i in 1:length(unique(gff[, 3]))) {
    message(unique(gff[, 3])[i])
  }
  
  message('\nCollecting signal...')
  message('(Skip genes with missing coordinates and signal in wiggle data)')

  # Create data frames to collect final data for all chrs
  plus_final <- data.frame()
  minus_final <- data.frame()
  # Keep track of total and non-skipped genes, to print info at the end
  number_genes <- 0
  number_skipped_genes <- 0
  
  for(i in 1:length(inputData)) {
    chrNum <- paste0('chr', chrom[i])
    message(paste0(chrNum, ':'))
    
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
      full_leng <- (end - start) + 1
      gene <- chromGff[j, 9]
      
      # Pull out signal
      sig_gene <- chromData[which(chromData[, 1] >= start & chromData[, 1] <= end), ]
      
      # Skip if there are discontinuities in the data (missing position:value pairs)
      if(nrow(sig_gene) != full_leng) next
      
      # Normalize to segment length of 1000
      sig_gene$position <- (sig_gene$position - start) + 1
      sig_gene$position <- sig_gene$position * (1000 / full_leng)
      
      # Genes of different sizes have different numbers of positions; small genes
      # (<1000bp) cannot produce signal values for all 1000 positions and will have gaps
      # This means that longer genes with more signal values per each position in the
      # sequence of 1000 positions will contribute more to the final output.
      # In order to avoid this, first build a Loess model with low span argument (alpha)
      # of the signal and then project it onto 1000 positions by using the model
      # to predict the signal
      model <- loess(signal ~ position, sig_gene,
                     control = loess.control(surface = "direct"),
                     span = loessSpan)
      signal <- predict(model, data.frame(position = seq(1, 1000, 1)), se = F)
      sig_gene <- data.frame(position = seq(1, 1000, 1), signal)
      
      # Make data frame for this gene
      all <- data.frame(chr = paste0('chr', chrom[i]), sig_gene, gene)

      # To collect all genes
      plus_sigs <- dplyr::bind_rows(plus_sigs, all)
      
      geneCount <- geneCount + 1
    }
    message(paste0('... + strand: ', geneCount, ' genes (skipped ', j - geneCount, ')'))
    
    # Keep track of total and non-skipped genes, to print info at the end
    number_genes <- number_genes + j
    number_skipped_genes <- number_skipped_genes + (j - geneCount)
    
    # To collect all chrs
    plus_final <- dplyr::bind_rows(plus_final, plus_sigs)
    
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
      full_leng <- (end - start) + 1
      gene <- chromGff[j, 9]
      
      # Pull out signal
      sig_gene <- chromData[which(chromData[, 1] >= start & chromData[, 1] <= end), ]
      
      # Skip if there are discontinuities in the data (missing position:value pairs)
      if(nrow(sig_gene) != full_leng) next
      
      # Normalize to segment length of 1000
      sig_gene$position <- (sig_gene$position - start) + 1
      sig_gene$position <- sig_gene$position * (1000 / full_leng)
      # Genes of different sizes have different numbers of positions; small genes
      # (<1000bp) cannot produce signal values for all 1000 positions and will have gaps
      # This means that longer genes with more signal values per each position in the
      # sequence of 1000 positions will contribute more to the final output.
      # In order to avoid this, first build a Loess model with low span argument (alpha)
      # of the signal and then project it onto 1000 positions by using the model
      # to predict the signal
      model <- loess(signal ~ position, sig_gene,
                     control = loess.control(surface = "direct"),
                     span = loessSpan)
      signal <- predict(model, data.frame(position = seq(1, 1000, 1)), se = F)
      sig_gene <- data.frame(position = seq(1, 1000, 1), signal)
      
      # Reverse the order of the position values
      sig_gene$position <- (1000 - sig_gene$position) + 1
      
      # Make data frame for this gene
      all <- data.frame(chr = paste0('chr', chrom[i]), sig_gene, gene)

      # To collect all genes
      minus_sigs <- dplyr::bind_rows(minus_sigs, all)
      
      geneCount <- geneCount + 1
    }
    # To collect all chrs
    minus_final <- dplyr::bind_rows(minus_final, minus_sigs)
    
    message(paste0('... - strand: ', geneCount, ' genes (skipped ', j - geneCount, ')\n'))
    # Keep track of total and non-skipped genes, to print info at the end
    number_genes <- number_genes + j
    number_skipped_genes <- number_skipped_genes + (j - geneCount)
  }
  
  # Merge '+' and '-' strand data
  mergedStrands <- dplyr::bind_rows(plus_final, minus_final)
  colnames(mergedStrands) <- c("chrom", "position", "signal", "gene")
  # Sort by gene and position
  mergedStrands <- mergedStrands[order(mergedStrands$gene, mergedStrands$position), ]
  
  message(paste0('Completed in ', round((proc.time()[3] - ptm[3]) / 60, 2), ' min.'))
  
  # Print info on total and non-skipped genes
  message('------')
  message(paste0('Skipped ', number_skipped_genes, ' of a total of ', number_genes,
             " genes (", round((number_skipped_genes * 100 / number_genes), 1),
             "%)."))
  message('------')
  
  if(saveFile) {
    message(paste0('Saving file...'))
    if(check_S288C) {
      write.table(mergedStrands, paste0(deparse(substitute(inputData)),
                                        "_S288C_metaORF.txt"), sep = "\t", quote = FALSE,
                                        row.names = FALSE)
      message('Done!')
    } else {
      write.table(mergedStrands, paste0(deparse(substitute(inputData)),
                                        "_SK1_metaORF.txt"), sep = "\t", quote = FALSE,
                                        row.names = FALSE)
      message('Done!')
    }
  } else {
    message('Done!')
    return(mergedStrands)
  }
}
