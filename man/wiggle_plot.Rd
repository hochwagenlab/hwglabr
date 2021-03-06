% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wiggle_plot.R
\name{wiggle_plot}
\alias{wiggle_plot}
\title{Wiggle data line plot}
\usage{
wiggle_plot(wiggleData, chr, genome, yMax, color = "grey50", protein,
  legendXcoord = -10, legendYcoord = yMax,
  legendAnnotation = deparse(substitute(wiggleData)), onScreen = TRUE)
}
\arguments{
\item{wiggleData}{A data frame of wiggle data with two columns: genome position and signal. No default.}

\item{chr}{A number representing the chromosome of \code{wiggleData}. No default.}

\item{genome}{A string representing the genome used for mapping. No default.}

\item{yMax}{Optional number to be used as the max Y scale value in the plots.
Particularly useful to plot two chromosomes on the same Y scale. No default.}

\item{color}{Optional R color. Defaults to \code{grey50}.}

\item{protein}{A string representing the ChIPped protein. No default.}

\item{legendXcoord}{A number representing the X coordinate to locate legend.
Defaults to minimum X (left-aligned).}

\item{legendYcoord}{A number representing the Y coordinate to locate legend.
Defaults to maximum Y (top-aligned).}

\item{legendAnnotation}{Optional string to be used as the legend. Defaults to name of object
passed to the function.}

\item{onScreen}{Boolean indicating plots should be returned to the screen (\code{onScreen = TRUE})
or written to .png files (\code{onScreen = FALSE}). Defaults to \code{TRUE}.}
}
\value{
A line plot, either on screen or as a png file (in the working directory).
}
\description{
This function allows you to plot wiggle data for a selected chromosome.
It takes as input either the raw wiggle data (for example one element of the 16 chromosome list)
or the output of \code{wiggle_smooth()} (R data frame with two columns: genome position and
smoothed signal).
}
\examples{
\dontrun{
wiggle_plot(WT[[1]], 1, genome = 'SK1', protein = 'Red1')

wiggle_plot(WT_chr3, 3, genome = 'SK1', yMax = 5, color = 'red', protein = 'Red1', onScreen = TRUE)

wiggle_plot(chrXVI, 16, genome = 'S288C', yMax = 5, color = 'black',
            protein = 'Rec8-HA', legendXcoord = 600, onScreen = FALSE)
}
}
