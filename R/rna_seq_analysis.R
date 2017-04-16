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
#' experiments without biological replicate libraries. While the calculation of CPM and
#' TPM has no requirement for replicates, you should obviously avoid doing DE analysis
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
#' @param batchNames Optional list of strings corresponding to the experimental batch for each
#' sample/library. Will be part of the design matrix for differential expression testing and
#' will define batch groups, if any. Defaults to \code{NULL} (all samples come from the same batch).
#' @param pairwiseDE Logical indicating whether to test differential expression (DE).
#' Defaults to \code{FALSE}.
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
#'   \item \strong{DE:} A \strong{D}ifferential \strong{E}xpression table comparing the first sample
#'   to all others. Includes log2 fold change (\code{logFC}), average log2 counts per million
#'   (\code{logCPM}), two-sided p-value (\code{PValue}) and false discovery rate (\code{FDR}).\cr
#'   \cr\strong{Note:} Fold change differences (\code{logFC}) between samples not directly ompared
#'   in the table can be obtained by subtracting their reported \code{logFC} (to the first sample).
#'   For example, if the first sample is \code{sample1} and we want the \code{logFC} between
#'   \code{sample2} and \code{sample3}, simply calculate the difference:\cr
#'   \cr\code{logFC.sample3_vs_sample2} = \code{logFC.sample3} - \code{logFC.sample2}
#' }
#' A \strong{M}ulti-\strong{D}imensional \strong{S}caling plot of all samples is also
#' saved in the output directory as a .pdf file.
#' @examples
#' \dontrun{
#' rna_seq_analysis(pathToFiles = list('AH119-2h_featureCounts.txt', 'AH119-3h_featureCounts.txt',
#'                                     'AH8104-2h_featureCounts.txt', 'AH8104-3h_featureCounts.txt'),
#'                  sampleNames = list('AH119_2h', 'AH119_3h', 'AH8104_2h', 'AH8104_3h'),
#'                  conditionNames = list('WT_2h', 'WT_3h', 'dot1_2h', 'dot1_3h'),
#'                  outputFilePrefix = 'dot1_noReplicates')
#'
#' rna_seq_analysis(pathToFiles = list('AH119-2h_featureCounts.txt', 'AH119-3h_featureCounts.txt',
#'                                     'AH8104-2h_featureCounts.txt', 'AH8104-3h_featureCounts.txt'),
#'                  sampleNames = list('AH119_2h', 'AH119_3h', 'AH8104_2h', 'AH8104_3h'),
#'                  conditionNames = list('WT_2h', 'WT_3h', 'dot1_2h', 'dot1_3h'),
#'                  pairwiseDE = TRUE, outputFilePrefix = 'dot1_noReplicates')
#'
#' rna_seq_analysis(pathToFiles = list('AH119-A_featureCounts.txt', 'AH119-B_featureCounts.txt',
#'                                     'AH8104-A_featureCounts.txt', 'AH8104-B_featureCounts.txt'),
#'                  sampleNames = list('AH119_A', 'AH119_B', 'AH8104_A', 'AH8104_B'),
#'                  conditionNames = list('WT', 'WT', 'dot1', 'dot1'),
#'                  batchNames = list('batch1', 'batch1', 'batch2', 'batch2')
#'                  pairwiseDE = TRUE)
#' }
#' @export

rna_seq_analysis <- function(pathToFiles, sampleNames, conditionNames, batchNames = NULL,
                             pairwiseDE = FALSE, outputFilePrefix){
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
  message('Created output directory "RNA-seq_analysis"')
  
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
  message('Created DGEList object:')
  print(y$samples[, 1:2])
  
  #----------------------------------------------------------------------------#
  #----------------------------- CPM and TPM ----------------------------------#
  # Filter out genes with low counts (low or no expression) based on CPMs,
  # in order to account for differences in library depth
  # Use a low threshold (cpm of 0.5)
  message('Filtering out features with counts below threshold (cpm < 0.5):')
  # If there are no replicates, keep all genes expressed in at least one sample
  # If there are replicates, keep all genes expressed in at least two samples
  # The code to filter rows in the edgeR manual no longer seems to work;
  # try it and fallback on an alternative in case it doesn't work
  try(
    if(length(unique(conditionNames)) == length(unique(sampleNames))){
      keep <- rowSums(edgeR::cpm(y) > 0.5) >= 1
      y_filt <- y[keep, , keep.lib.sizes = FALSE]
    } else {
      keep <- rowSums(edgeR::cpm(y) > 0.5) >= 2
      y_filt <- y[keep, , keep.lib.sizes = FALSE]
    },
  silent = T)
 
  if(!exists('y_filt')){
    y_filt <- y
    y_filt$counts <- y$counts[keep, ]
    y_filt$genes <- y$genes[keep, ]
  }
  
  # Calculate updated lib sizes (differences should be minimal):
  y_filt$samples$lib.size <- colSums(y_filt$counts)
   
  message(nrow(y_filt$counts), ' features kept of an original total of ',
          nrow(y$counts))
  dropped <- round((nrow(y$counts) - nrow(y_filt$counts)) / nrow(y$counts) * 100, 1)
  message(' (', dropped, '% filtered out)')
  
  # Calculate normalization factors to scale the raw library sizes
  y_filt <- edgeR::calcNormFactors(y_filt)
  message('Calculated normalization factors using trimmed mean of M-values (TMM) method.')
  y_filt$samples
  
  # Save MDS plot
  if(missing(outputFilePrefix)){
    pdf(paste0("RNA-seq_analysis/", "MDS_plot.pdf"))
  } else{
    pdf(paste0("RNA-seq_analysis/", outputFilePrefix, "_MDS_plot.pdf"))
  }
  limma::plotMDS(y_filt)
  dev.off()
  message('Plotted multi-dimensional scaling (MDS) and saved to .pdf file.')
  
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
  message('Calculated CPM and saved to file.')
  
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
  message('Calculated TPM and saved to file.')
  
  #----------------------------------------------------------------------------#
  #----------------------------- DE analysis ----------------------------------#
  # Run if not missing
  if(pairwiseDE){
    de <- de_analysis(DGEListObject=y_filt, sampleNames, conditionNames,
                      batchNames, outputFilePrefix)
    
    if(missing(outputFilePrefix)){
      write.csv(de$table, paste0("RNA-seq_analysis/", "de.csv"), row.names = T)
    } else {
      write.csv(de$table, paste0("RNA-seq_analysis/", outputFilePrefix, "_de.csv"),
                row.names = T)
    }
    message('Calculated DE and saved to file.')
  }
  
  message("... ... ...")
  message("(Output files are in directory \"RNA-seq_analysis\").")
  message("Completed in ", round((proc.time()[3] - ptm[3]), 1), "sec.")
}



### Helper function to perform differential expression analyses
# Note: The experimental design is parametrized with a one-way layout.
# Must not be used if this is not appropriate for the analyzed experiment
de_analysis <- function(DGEListObject, sampleNames, conditionNames,
                        batchNames, outputFilePrefix){
  message('Running DE analyses:')
  # Design matrix
  group <- factor(unlist(conditionNames))
  if(!is.null(batchNames)){
    batch <- factor(unlist(batchNames))
    design <- model.matrix(~batch+group)
    colnames(design)[2:ncol(design)] <- c(tail(levels(batch), -1),
                                          tail(levels(group), -1))
  } else {
    design <- model.matrix(~group)
    colnames(design)[2:ncol(design)] <- tail(levels(group), -1)
  }
  
  ### Test DE:
  # While the likelihood ratio test is a more obvious choice for inferences with
  # GLMs, the QL F-test is preferred as it provides more robust and reliable
  # error rate control when the number of replicates is small.
  # glmQLFit() can only be used when there are replicates, however.
  # In the absence of replicates, use glmFit() followed by glmLRT() (using the
  # typical value for BCV - square root-dispersion).
  
  # Estimate dispersion and fit model
  if(length(unique(conditionNames)) == length(unique(sampleNames))){
    message('Performing DE analysis (without replicates!)')
    bcv <- 0.1
    fit <- edgeR::glmFit(DGEListObject, design, dispersion = bcv^2)
    
    startCondition <- ifelse(!is.null(batchNames),
                             length(levels(batch)) + 1, 2)
    de_test <- edgeR::glmLRT(fit, coef=startCondition:ncol(design))
  } else {
    message('\nPerforming DE analysis with replicates:\n')
    y_filt <- edgeR::estimateDisp(DGEListObject, design, robust=TRUE)
    message('Estimated common dispersion: ', y_filt$common.dispersion)
    fit <- edgeR::glmQLFit(y_filt, design)
    startCondition <- ifelse(!is.null(batchNames),
                             length(levels(batch)) + 1, 2)
    de_test <- edgeR::glmQLFTest(fit, coef=startCondition:ncol(design))
  }
  
  # Get all genes (n=Inf)
  de_test <- edgeR::topTags(de_test, n=Inf)
  return(de_test)
}
