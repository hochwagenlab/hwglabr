#' Wrapper function to run edgeR analysis of summarized reads from RNA-seq experiment
#'
#' This function allows you to run a statistical analysis of an RNA-seq experiment
#' using Bioconductor package
#' \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR}.
#' Inputs are tables of summarized reads generated using
#' \href{http://bioinf.wehi.edu.au/featureCounts/}{featureCounts} (part of the lab's
#' RNA-seq analysis pipeline running on NYU's HPC) for each sample in the experiment.
#' The presence or absence of biological replicates is automatically inferred from the
#' \code{conditionNames} argument: conditions aggregate samples as replicates. If
#' you provide as many different \code{conditionNames} as there are \code{sampleNames}
#' each sample is taken as a single condition (no replicates).
#' The output includes tables of CPM, TPM and, if included, differential expression
#' (DE) analysis for selected pairs of samples (see "Value" section below for more
#' details).\cr
#' \cr
#' \strong{Running without replicates:} This function allows you to run the analysis on
#' experiments without replicate libraries. While the calculation of CPM and TPM
#' has no requirement for replicates, you should obviously avoid doing DE analysis
#' in the absence of replicates. The statistical model used by
#' \href{https://bioconductor.org/packages/release/bioc/html/edgeR.html}{edgeR}
#' to perform DE analysis relies on biological replicates to estimate biological
#' variability. While there is no satisfactory alternative to actual biological replication,
#' you can still run a preliminary DE analysis without replicates. The approach
#' taken here relies on providing an estimate of biological variation that is reasonable
#' for RNA-seq experiments using \emph{S. cerevisiae}. This is a better alternative
#' to simply assuming that biological variability is absent. Typical values for
#' the common BCV (square root-dispersion) for datasets arising from well-controlled
#' experiments are 0.1 for data on genetically identical model organisms, so that
#' value is used here.
#' @param pathToFiles A list of strings corresponding to the full path to the featureCounts
#' output files for all samples in the experiment. No default.
#' @param sampleNames A list of strings corresponding to the names of all samples
#' in the experiment in the same order as the files in \code{pathToFiles}. No default.
#' @param conditionNames A list of strings corresponding to the experimental groups/conditions
#' for each sample/library. Will be the input to edgeR's \code{\link[edgeR]{DGEList}}
#' constructor function and will define replicate groups, if any. No default.
#' @param pairwiseDE Optional list of two-element string vectors corresponding to
#' the names of pairs of samples to test differential expression (DE) for. Must
#' match strings in \code{conditionNames}. No default.
#' @param outputFilePrefix Optional string to be added as prefix to output file names.
#' No default.
#' @return The output includes several tables saved as .csv files in a directory
#' named "RNA-seq_analysis" written to the working directory:
#' \enumerate{
#'   \item \strong{CPM:} Compositional bias-normalized \strong{C}ounts \strong{P}er
#'   \strong{M}illion (output of edgeR's \code{\link[edgeR]{cpm}} function). Useful
#'   for intersample feature expression comparison (not feature length-normalized).
#'   \item \strong{TPM:} Compositional bias and feature length-normalized \strong{T}ranscripts
#'   \strong{P}er \strong{M}illion. Useful for within-sample feature expression comparison.
#'   \item \strong{DE:} A \strong{D}ifferential \strong{E}xpression table for each
#'   requested pair of samples (output of edgeR's \code{\link[edgeR]{exactTest}} function),
#'   including log2 fold change (\code{logFC}), average log2 counts per million
#'   (\code{logCPM}), and two-sided p-value (\code{PValue}).
#' }
#' A \strong{M}ulti-\strong{D}imensional \strong{S}caling plot of all samples is also
#' saved in the output directory as a .pdf file.
#' @examples
#' \dontrun{
#' rna_seq_analysis(pathToFiles = list('AH119-2h_featureCounts.txt', 'AH119-3h_featureCounts.txt',
#'                                     'AH8104-2h_featureCounts.txt', 'AH8104-3h_featureCounts.txt'),
#'                  sampleNames = list('AH119_2h', 'AH119_3h', 'AH8104_2h', 'AH8104_3h'),
#'                  conditionNames = list('WT_2h', 'WT_3h', 'dot1_2h', 'dot1_3h'),
#'                  pairwiseDE = list(c('WT_2h', 'dot1_2h'), c('WT_3h', 'dot1_3h'),
#'                                    c('WT_2h', 'dot1_2h'), c('WT_3h', 'dot1_3h'),
#'                                    c('WT_2h', 'WT_3h'), c('dot1_2h', 'dot1_3h')),
#'                  outputFilePrefix = 'dot1_noReplicates')
#'
#'  rna_seq_analysis(pathToFiles = list('AH119-2h_featureCounts.txt', 'AH119-3h_featureCounts.txt',
#'                                     'AH8104-2h_featureCounts.txt', 'AH8104-3h_featureCounts.txt'),
#'                  sampleNames = list('AH119_2h', 'AH119_3h', 'AH8104_2h', 'AH8104_3h'),
#'                  conditionNames = list('WT_2h', 'WT_3h', 'dot1_2h', 'dot1_3h'),
#'                  outputFilePrefix = 'dot1_noReplicates')
#'
#' rna_seq_analysis(pathToFiles = list('AH119-A_featureCounts.txt', 'AH119-B_featureCounts.txt',
#'                                     'AH8104-A_featureCounts.txt', 'AH8104-B_featureCounts.txt'),
#'                  sampleNames = list('AH119_A', 'AH119_B', 'AH8104_A', 'AH8104_B'),
#'                  conditionNames = list('WT', 'WT', 'dot1', 'dot1'),
#'                  pairwiseDE = list(c('WT', 'dot1')))
#' }
#' @export

rna_seq_analysis <- function(pathToFiles, sampleNames, conditionNames,
                             pairwiseDE, outputFilePrefix){
  ptm  <- proc.time()
  
  #----------------------------------------------------------------------------#
  #-------------------------- Preliminary checks ------------------------------#
  if (file.exists('RNA-seq_analysis')) {
    stop('ERROR: A folder called "RNA-seq_analysis" already exists in the current working directory.\n',
         'Please remove it and repeat function call.', call. = FALSE)
  }
  
  if (!all(unlist(lapply(pathToFiles, file.exists)), na.rm = FALSE)) {
    stop('ERROR: Could not open one or more featureCount files.',
         'Please check the provided paths to the files.', call. = FALSE)
  }
  
  if (!missing(pairwiseDE)) {
    if (!all(unlist(lapply(pairwiseDE, is.element, conditionNames)), na.rm = FALSE)) {
      stop('ERROR: strings for "pairwiseDE" must match "conditionNames".', call. = FALSE)
    }
  }
  
  # Check for dplyr and edgeR (and load edgeR)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("R package 'dplyr' is required. Please install it.\n",
         "install.packages('dplyr')", call. = FALSE)
  }
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Bioconductor package 'edgeR' is required. Please install it:\n", 
         "## try http:// if https:// URLs are not supported",
         "source('https://bioconductor.org/biocLite.R')",
         "biocLite('edgeR')", call. = FALSE)
  }
  
  # Create directory in current working directory to save output
  dir.create('RNA-seq_analysis')
  message('Created output directory "RNA-seq_analysis"\n')
  
  #----------------------------------------------------------------------------#
  #----------------------- Generate DGEList object ----------------------------#
  # Load count data and create table of all samples
  # Load files
  counts <- lapply(pathToFiles, read.table, header = T)
  # Rename feature counts column with supplied sample names
  for(i in 1:length(counts)){
    colnames(counts[[i]])[7] <- sampleNames[i]
  }
  # Reduce list to data frame
  counts <- Reduce(dplyr::inner_join, counts)
  # Create DGEList object
  y <- edgeR::DGEList(counts = counts[, unlist(sampleNames)],
                      group = unlist(conditionNames),
                      genes = data.frame(Gene = counts$Geneid,
                                         Length = counts$Length))
  message('Created DGEList object:\n')
  print(y$samples[, 1:2])
  
  #----------------------------------------------------------------------------#
  #----------------------------- CPM and TPM ----------------------------------#
  # Filter out genes with low counts (low or no expression) based on CPMs,
  # in order to account for differences in library depth
  # Use a low threshold (cpm of 0.5)
  message('Filtering out features with counts below threshold (cpm < 0.5):\n')
  # If there are no replicates, keep all genes expressed in at least one sample
  # If there are replicates, keep all genes expressed in at least two samples
  if(length(unique(conditionNames)) == length(unique(sampleNames))){
    keep <- rowSums(edgeR::cpm(y) > 0.5) >= 1
    y_filt <- y[keep, , keep.lib.sizes = FALSE]
  } else {
    keep <- rowSums(edgeR::cpm(y) > 0.5) >= 2
    y_filt <- y[keep, , keep.lib.sizes = FALSE]
  }
  
  message(nrow(y_filt$counts), 'features kept of an original total of', nrow(y$counts))
  dropped <- round((nrow(y$counts) - nrow(y_filt$counts)) / nrow(y$counts) * 100, 1)
  message(' (', dropped, '% filtered out).\n', sep = '')
  
  # Calculate normalization factors to scale the raw library sizes
  y_filt <- calcNormFactors(y_filt)
  message('Calculated normalization factors using trimmed mean of M-values (TMM) method:\n')
  message(y_filt$samples)
  
  # Save MDS plot
  if(missing(outputFilePrefix)){
    pdf(paste0("RNA-seq_analysis/", "MDS_plot.pdf"))
  } else{
    pdf(paste0("RNA-seq_analysis/", outputFilePrefix, "_MDS_plot.pdf"))
  }
  limma::plotMDS(y_filt)
  dev.off()
  message('Plotted multi-dimensional scaling (MDS) and saved to .pdf file.\n')
  
  # Calculate CPM and write to file
  cpm_edgeR <- edgeR::cpm(y_filt)
  rownames(cpm_edgeR) <- y_filt$genes$Gene
  
  # Write to file
  if(missing(outputFilePrefix)){
    write.csv(cpm_edgeR, paste0("RNA-seq_analysis/", "edgeR_cpm.csv"), row.names = T)
  } else {
    write.csv(cpm_edgeR, paste0("RNA-seq_analysis/", outputFilePrefix, "_edgeR_cpm.csv"),
              row.names = T)
  }
  message('Calculated CPM and saved to file.\n')
  
  # Calculate TPM
  # Write tpm function
  calcTPM <- function(inputDGEList, gene.length) {
    x <- as.matrix(inputDGEList)
    len.norm.lib.size <- colSums(x / gene.length)
    tpm_pre_len_norm <- t(t(x) / len.norm.lib.size) * 1e06
    rownames(tpm_pre_len_norm) <- inputDGEList$genes$Gene
    return(tpm_pre_len_norm / gene.length)
  }
  
  tpm <- calcTPM(y_filt, gene.length = y_filt$genes$Length)
  
  
  ### Write TPMs to file
  # Prep table
  rownames(tpm) <- y_filt$genes$Gene
  
  # Write to file
  if(missing(outputFilePrefix)){
    write.csv(tpm, paste0("RNA-seq_analysis/", "tpm.csv"), row.names = T)
  } else {
    write.csv(tpm, paste0("RNA-seq_analysis/", outputFilePrefix, "_tpm.csv"), row.names = T)
  }
  message('Calculated TPM and saved to file.\n')
  
  #----------------------------------------------------------------------------#
  #----------------------------- DE analysis ----------------------------------#
  # Run if not missing
  if(!missing(pairwiseDE)){
    if(length(unique(conditionNames)) == length(unique(sampleNames))){
      message('\nPerforming DE analysis (without replicates!):\n')
      bcv <- 0.1
    } else {
      message('\nPerforming DE analysis with replicates:\n')
      y_filt <- estimateDisp(y_filt)
    }
    
    for(i in 1:length(pairwiseDE)){
      if(length(unique(conditionNames)) == length(unique(sampleNames))){
        et <- edgeR::exactTest(y_filt, pair = pairwiseDE[[i]], dispersion = bcv^2)
      } else {
        et <- edgeR::exactTest(y_filt, pair = pairwiseDE[[i]])
      }
      rownames(et$table) <- y_filt$genes$Gene

      if(missing(outputFilePrefix)){
        write.csv(et$table, paste0("RNA-seq_analysis/DE_", et$comparison[1], "-",
                                   et$comparison[2], ".csv"), row.names = T)
      } else {
        write.csv(et$table, paste0("RNA-seq_analysis/", outputFilePrefix, "_DE_",
                                   et$comparison[1], "-", et$comparison[2], ".csv"),
                  row.names = T)
      }
      message(et$comparison[1], "vs", et$comparison[2],
          "saved to file.\n")
    }
  }
  
  message("... ... ...\n(Output files are in directory \"RNA-seq_analysis\").\n")
  message("Completed in", round((proc.time()[3] - ptm[3]), 1), "sec.")
}
