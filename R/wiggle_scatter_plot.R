#' To make a scatter plot comparison of two wiggle files
#'
#' This function is designed to work with function \code{\link{wiggle_scatter}} to
#' allow for an overall generalized comparison of two data sets.
#' Great for comparing replicates, checking controls, etc.
#' The two functions are used as follows:
#' \enumerate{
#'   \item \code{\link{wiggle_scatter}} will extract the necessary information using
#'   \code{\link{wiggle_compress}} for all chromosomes
#'   \item \code{wiggle_scatter_plot()} will make a simple .eps file of the results.
#' }
#' Written by Tovah Markowitz
#' @param wiggleData1 Parameter of \code{\link{wiggle_scatter}}. A list of 16 chr wiggle data 
#' (output of \code{\link{readall_tab}}. No default.
#' @param wiggleData2 Parameter of \code{\link{wiggle_scatter}}. A second set of wiggle data 
#' (output of \code{readall_tab()}). No default.
#' @param window Parameter of \code{\link{wiggle_scatter}}. Size of window to be compressed.
#' Default is 5000 bp.
#' @param scatter Parameter for \code{\link{wiggle_scatter}}. Input should be output of 
#' \code{wiggle_scatter()}.
#' @return Output of \code{\link{wiggle_scatter}} is a large (dplyr) data frame with
#' compressed data from all chromosomes. See \code{\link{wiggle_compress}} for more 
#' information. \cr
#' Output of \code{wiggle_scatter_plot()} is an .eps file with file name and axis
#' names determined by columns of input data frame. The function is written so that
#' the x-axis and the y-axis have the same range.
#' @examples
#' a <- wiggle_scatter(red1, hop1, 4000)
#' 
#' wiggle_scatter_plot(a)
#' @export

wiggle_scatter_plot <- function( scatter ) {
  # example scatter plot function, made as eps for use in illustrator
  LIM <- range( c( scatter[,2], scatter[,3] ), na.rm = T)
  
  setEPS()
  postscript( paste0( names(scatter)[2], "vs", names(scatter)[3], "scatterplot.eps") )
  plot(unlist(scatter[,2]), unlist(scatter[,3]), xlim = LIM, ylim = LIM,
       pch = 20, xlab = names(scatter)[2], ylab = names(scatter)[3] )
  dev.off()
}