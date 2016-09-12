#' To make a scatter plot comparison of two wiggle files
#'
#' This function is designed to work with function \code{\link{wiggle_scatter}} to
#' allow for an overall generalized comparison of two data sets.
#' Great for comparing replicates, checking controls, etc.
#' The two functions are used as follows:
#' \enumerate{
#'   \item \code{\link{wiggle_scatter}} will extract the necessary information using
#'   \code{\link{wiggle_compress}} for all chromosomes
#'   \item \code{wiggle_scatter_plot} will make a simple .eps file of the results.
#' }
#' Written by Tovah Markowitz
#' @param scatter Input should be output of \code{\link{wiggle_scatter}}.
#' @param method Method used to calculate correlation value on plot (same as for
#' \code{\link{cor}} function). Defaults to "spearman".
#' @return Output is an .eps file with file name and axis
#' names determined by columns of input data frame. The function is written so that
#' the x-axis and the y-axis have the same range.
#' @examples
#' \dontrun{
#' a <- wiggle_scatter(red1, hop1, 4000)
#' 
#' wiggle_scatter_plot(a)
#' }
#' @export

wiggle_scatter_plot <- function( scatter, method = "spearman" ) {
  # example scatter plot function, made as eps for use in illustrator
  LIM <- range( c( scatter[,2], scatter[,3] ), na.rm = T)
  
  corr <- cor(unlist(scatter[,2]),unlist(scatter[,3]),method=method,use="complete.obs")
  
  setEPS()
  postscript( paste0( names(scatter)[2], "vs", names(scatter)[3], "scatterplot.eps") )
  plot(unlist(scatter[,2]), unlist(scatter[,3]), xlim = LIM, ylim = LIM,
       pch = 20, xlab = names(scatter)[2], ylab = names(scatter)[3] )
  text(LIM[1]+0.5,LIM[2],bquote(rho==.(corr)))
  text(LIM[1]+0.5,LIM[2]-0.5,bquote(R^2==.(corr^2)))
  dev.off()
}