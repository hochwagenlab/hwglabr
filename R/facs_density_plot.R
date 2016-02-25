#' FACS density plot
#'
#' This function allows you to plot a density plot across time points of FACS data.
#' Written by Tovah Markowitz.
#' @param yeastLine A number present in the file names corresponding to the yeast
#' strain number. No default.
#' @param gates Optional gating limits as \code{c(gateLower, gateUpper)}. Defaults to no gating.
#' @param type A string representing the image file type, acepting either \code{jpg} or \code{pdf}.
#' Defaults to \code{jpg}.
#' @section Details:
#' The function is designed to search all .FCS files (named as "Cell YeastLine_TimePoint.fcs")
#' for those with a particular yeast line and print a quick view of ungated FACS
#' results to the screen. It also allows for the option of gating the data before
#' writing the results to an image file in the working directory.
#' 
#' The function requires two R packages: "flowCore" and "flowViz".
#' You can install the packages by running:
#' 
#' \code{source("http://bioconductor.org/biocLite.R")}
#' \code{biocLite("flowCore")}
#' \code{biocLite("flowViz")}
#' 
#' (If the packages are missing you will see an error message with this information).
#' @section Note:
#' If the number of the yeast line is found anywhere else in the name of any of
#' the files in the working directory, this function WILL NOT work.
#' @return A density plot of all time points with the same yeast strain,
#' either on screen or as a file (in the working directory).
#' @examples
#' facs_density_plot(119)
#' 
#' facs_density_plot(119, c(2000000, 7500000), 'pdf')
#' @export



facs_density_plot <- function(yeastLine, gates, type) {
  
  # This piece allows the gating to be optional
  # along with the other 'if statement' later in the code
  if (missing(gates)) {
    gates = vector()
  }
  if (missing(type)) {
  	type = "jpg"
  }
  # Open the required packages
  if (!requireNamespace("flowCore", quietly = TRUE)) {
    stop("R package 'flowCore' needed for this function to work. Please install it:\n",
         "source('http://bioconductor.org/biocLite.R')\n", "biocLite('flowCore')",
         call. = FALSE)
  }
  if (!requireNamespace("flowViz", quietly = TRUE)) {
    stop("R package 'flowViz' needed for this function to work. Please install it:\n",
         "source('http://bioconductor.org/biocLite.R')\n", "biocLite('flowViz')",
         call. = FALSE)
  }
  
  library("flowCore")
  library("flowViz")
  # Extract the files with the correct yeastLine. 
  # Note: if the number of the yeast line is found anywhere else in the name 
  # of any of the files in the working directory this function WILL NOT work. 
  # Yeast lines to watch: 1-31. 
  
  files <- list.files()
  yfiles <- files[grep(yeastLine, files)]
  fileNames <- yfiles[grep("fcs", yfiles)]
  numFiles <- length(fileNames)
  
  # For each file found:
  samples <- vector("list")
  for (z in 1:numFiles) {
    # Determine the time point of the file. The internal string split accesses 
    # everything before the period [[1]][1]. The external string split accesses 
    # the third subset of the file name [[1]][3]. This assumes that the format of 
    # the file name is: date_yeastLine_timePoint.fcs
    timePoint <- strsplit(strsplit(fileNames[z], "\\.") [[1]][1], "_")[[1]][2]
    
    # Read the FCS file
    samples[[timePoint]] <- flowCore::read.FCS(fileNames[z])
  }
  
  # Create a FlowSet from the samples across time points and create a quick plot
  # of the FITC-A results. To change the features of the density plot, read more
  # about FlowCore, FlowViz, Lattice, and Trellis
  # Note: densityplot is a Lattice plot, not a standard R plot.
  
  experiment <- as(samples, "flowSet")
  print(densityplot(~`FL1-A`, experiment, overlap = 0.8))
  
  if (length(gates) != 0) {
    
    # To use rectangular gating to filter the data. Information on filterIds is
    # limited. Note: some smoothing also occurs with this function
    rectGate <- rectangleGate(filterId="nonDebris","FL1-A" = gates)
    filtered <- filter(experiment,rectGate)
    
    if (type == "jpg") {
      fileName <- as.character(paste(yeastLine, ".jpg", sep = ""))
      # This is the standard way to change Lattice features. Any features of the
      # plot that cannot be changed using basic XY-plot features should be changed
      # in this way; otherwise, all future Lattice plots may continue to have these
      # features.
      
      # Step 1: Open a device in which to save the Lattice
      trellis.device(jpeg, file = fileName)
      } else if (type == "pdf") {
        fileName <- as.character(paste(yeastLine, ".pdf", sep = ""))
        trellis.device(pdf, file = fileName)
      }
    
    # Step 2: Create the figure, but do not print. The fitc.plot call is essential
    # here.
    
    Xlow <- gates[1]-1000000
    Xhigh <- gates[2]+1000000
    fitc.plot <- densityplot(~`FL1-A`, Subset(experiment, filtered), overlap = 0.8,
                           main = list(as.character(yeastLine), fontsize = 20),
                           xlim = c(Xlow, Xhigh),
                           par.strip.text = list(cex = 1.5),
                           scales = list(cex = 2))
    
    # Step 3: Get the default settings for the feature of interest
    superpose.settings <- trellis.par.get("superpose.polygon")
    # Step 4: Replace the default settings with new settings for this device only.
    # In this case, we are making all of the curves light grey.
    superpose.settings$col <- "grey"
    trellis.par.set("superpose.polygon", superpose.settings)
    
    # Step 5: Print the figure for the first time directly to the JPEG file
    print(fitc.plot)
    
    # Close the device to prevent overlays. This will also cause the Lattice
    # features to revert to default.
    dev.off()
  }
}
