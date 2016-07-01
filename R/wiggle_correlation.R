#' Correlation between two ChIP-seq experiments
#' 
#' This function allows you to calculate correlations between wiggle files
#' for ChIP-seq, both on an individual chromosome basis and genome-wide. It
#' takes as an input two lists of 16 chromosomes (output of \code{\link{readall_tab}}).
#' @param Set1 A list of 16 chr wiggle data (output of \code{\link{readall_tab}})
#' @param Set2 A list of 16 chr wiggle data (output of \code{\link{readall_tab}})
#' @param method A character string indicating which correlation coefficient
#' is to be computed (\code{method} argument to R function \code{cor()}).
#' Accepts one of \code{pearson} (default), \code{kendall}, or \code{spearman}.
#' @return A vector with 17 values with names indicating the chromosome associated
#' with each individual correlation value.
#' @examples
#' wiggle_correlation(Red1WT_rep1, Red1WT_rep2)
#' 
#' wiggle_correlation(Red1WT_rep1, Red1WT_rep2, method = "spearman")
#' @export

# By: Tovah Markowitz
# Date: 5/18/16

wiggle_correlation <- function(Set1, Set2, method = 'pearson') {
	ptm <- proc.time()
	
	if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("R package 'dplyr' needed for this function to work. Please install it.\n", 
            "install.packages('dplyr')", call. = FALSE)
    }
    #library(dplyr)

	# read in wiggle files to compare
  #	Set1 <- readall_tab(file1)
  #	Set2 <- readall_tab(file2)
	
	# complete correlation for each chromosome and entire genome
	correlation <- vector(length=16)
	b <- data.frame("position"=numeric(), "set1"=numeric(), "set2"=numeric())
	for (f in 1:16) {
		names(Set1[[f]]) <- c( "position", "set1")
		names(Set2[[f]]) <- c( "position", "set2")
		a <- dplyr::inner_join( Set1[[f]], Set2[[f]] )
		b <- dplyr::bind_rows(b, a)
		correlation[f] <- cor( a$set1, a$set2, method = method )
	}
	correlation[17] <- cor( b$set1, b$set2, method = method )
	
	# to inform which chromosome matches which value
	tmp <- unlist( strsplit( unlist( strsplit(names(Set1), split=".wig") ) ,split="_") )
	names(correlation) <- c(tmp[ grepl("chr",tmp) ], "all")
	
	cat(paste0("Completed in ", round((proc.time()[3] - ptm[3])/60, 
        2), " min.\n"))
        
    return(correlation)
}