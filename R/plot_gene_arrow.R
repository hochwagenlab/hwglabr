#' Plot genes on an R wiggle plot
#'
#' This function allows you to you to add arrows representing genes to an R wiggle plot. \cr
#' Written by Tovah Markowitz.
#' 
#' @param geneEnd A numeric position on the x-axis of the 3' end of the gene. No default.
#' @param geneLength The length of the gene. No default.
#' @param orientation A value of 1 for genes on the Watson strand or -1 for genes on the Crick strand.
#' No default.
#' @param minima A minimal gene length before the arrow starts turning inside out.
#' Dependent on the range of the x-axis. Default is 400 bp.
#' @param color Define the color of the gene. Default is black.
#' @param height To adjust to the height of the gene bar. Default is 0.2.
#' @param yPos To adjust the position of the gene on the y-axis. Default is -1. 
#' @examples
#' plot_gene_arrow(Red1_end, Red1_end-Red1_start, 1)
#' plot_gene_arrow(672823, 672823-670340, 1)
#' @export

plot_gene_arrow <- function(geneEnd, geneLength, orientation, minima = 400,
                            color = "black", height = 0.2, yPos = -1) {
	if (geneLength > minima) {
			xc <- c(0, minima, geneLength, geneLength, minima, 0)
	} else {
		xc <- c(0, geneLength, geneLength, geneLength, geneLength, 0)
	}
	yc <- c(0, height/2, height/2, -(height/2),-(height/2), 0)
	polygon(geneEnd-orientation*(xc), (yPos)+(yc), col=color)
}