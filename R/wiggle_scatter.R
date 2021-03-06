#' To make a scatter plot comparison of two wiggle files
#'
#' This function is designed to work with function \code{\link{wiggle_scatter_plot}} to
#' allow for an overall generalized comparison of two data sets.
#' Great for comparing replicates, checking controls, etc.
#' The two functions are used as follows:
#' \enumerate{
#'   \item \code{wiggle_scatter} will extract the necessary information using
#'   \code{\link{wiggle_compress}} for all chromosomes
#'   \item \code{\link{wiggle_scatter_plot}} will make a simple .eps file of the results.
#' }
#' Written by Tovah Markowitz
#' @param wiggleData1 A list of 16 chr wiggle data (output of \code{\link{readall_tab}}). No default.
#' @param wiggleData2 A second set of wiggle data (output of \code{\link{readall_tab}}). No default.
#' @param window Parameter of \code{\link{wiggle_scatter}}. Size of window to be compressed.
#' Default is 5000 bp.
#' @return Output is a large (dplyr) data frame with
#' compressed data from all chromosomes. See \code{\link{wiggle_compress}} for more 
#' information.
#' @examples
#' \dontrun{
#' a <- wiggle_scatter(red1, hop1, 4000)
#' 
#' wiggle_scatter_plot(a)
#' }
#' @export

wiggle_scatter <- function( wiggleData1, wiggleData2, window = 5000 ) {
  
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
