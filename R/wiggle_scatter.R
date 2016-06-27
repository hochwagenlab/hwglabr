#' To make a scatter plot comparison of two wiggle files
#'
#' This function is designed to work with function \code{wiggle_scatter_plot()} to
#' allow for an overall generalized comparison of two data sets.
#' Great for comparing replicates, checking controls, etc.
#' The two functions are used as follows:
#' \enumerate{
#'   \item \code{wiggle_scatter()} will extract the necessary information using
#'   \code{wiggle_compress()} for all chromosomes
#'   \item \code{wiggle_scatter_plot()} will make a simple .eps file of the results.} \cr
#' Written by Tovah Markowitz
#' @param wiggleData1 Parameter of wiggle_scatter. A list of 16 chr wiggle data 
#' (output of readall_tab). No default.
#' @param wiggleData2 Parameter of wiggle_scatter. A second set of wiggle data 
#' (output of readall_tab). No default.
#' @param window Parameter of wiggle_scatter. Size of window to be compressed.
#' Default is 5000 bp.
#' @param scatter Parameter for plot_scatter. Input should be output of 
#' wiggle_scatter.
#' @return Output of wiggle_scatter is a large (dplyr) data frame with
#' compressed data from all chromosomes. See \code{wiggle_compress()} for more 
#' information. \cr
#' Output of plot_scatter is an .eps file with file name and axis names
#' determined by columns of input data frame. The function is written so that
#' the x-axis and the y-axis have the same range.
#' @examples
#' a <- wiggle_scatter( red1, hop1, 4000 )
#' plot_scatter( a )
#' @export


wiggle_scatter <- function( wiggleData1, wiggleData2, window=5000 ) {

    # prepare the two files for wiggle_compress
    for (j in 1:16) {
    	fileN1 <- deparse( substitute( wiggleData1 ) )
    	f1 <- paste0( fileN1, "[[", j, "]]" )
    	fileN2 <- deparse( substitute( wiggleData2 ) )
    	f2 <- paste0( fileN2, "[[", j, "]]" )
    	# run wiggle_compress for each chromosome and combine all data into a single data.frame
    	
    	if (j == 1) {
    	  data <- wiggle_compress( c(f1, f2), window )
    	  names(data) <- c("position", fileN1, fileN2 )
    	} else {
    	  tmp <- wiggle_compress( c(f1, f2), window )
	      names(tmp) <- c("position", fileN1, fileN2 )
	      data <- dplyr::bind_rows( data, tmp )
    	}
    }
    return(data)
}