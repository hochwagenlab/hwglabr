#' Compress wiggle data
#'
#' This function will compress wiggle data by averaging values (including position) across
#' non-overlapping windows. It will take any number of wiggle files, but can
#' only handle one chromosome or chromosomal region per function call. This
#' function assumes that all data is from the chromosome and ignores all 
#' chromosomal information. This function ignores all positions without reads mapping
#' which may lead to slight variations in position depending on the wiggle data being used.
#' Note: calls function compress. \cr
#' Written by Tovah Markowitz.
#' @param inputWiggles Will accept two types of wiggle information:
#' \enumerate{
#'   \item A single wiggle region/chromosome
#'   \item Multiple wiggle files to be directly compared/compressed. These must be
#'   listed as strings in the following format:
#'   \code{c("wiggleData1[[chr]]", "wiggleData2[[chr]]", ...)}.
#' }
#' @param window Size of window to be averaged across. Suggested size for
#' chromosome-wide plots for presenting and publication: 200 bp. No default.
#' @return A R (dplyr) data frame with averaged position and value for each 
#' wiggle. Column names are defined by input strings/variables.
#' @examples
#' \dontrun{
#' wiggle_compress( red1[[6]], 200 )
#'
#' wiggle_compress( c( "red1[[6]]", "hop1[[6]]", "smc4[[6]]" ), 200 )
#' }
#' @export
 
wiggle_compress <- function( inputWiggles, window ) {
	message( "Compressing ", length(inputWiggles), " wiggle files...")
	ptm <- proc.time()

	if (!requireNamespace("dplyr", quietly = TRUE)) {
           stop("R package 'dplyr' needed for this function to work. Please install it.\n", 
            "install.packages('dplyr')", call. = FALSE)
    	}

	# if only one wiggle file
	if ( length(inputWiggles) == 1 ) {
	   # allows both character or variable if a single wiggle file
	   if ( is.character(inputWiggles) ) {
	      	data <- eval( parse( text = inputWiggles ) )
	      	names(data) <- c( "position", inputWiggles )
	   } else { stop("Incorrect input format. Check documentation for details.") }
           
	} else if ( is.data.frame(inputWiggles) ) {
           if ( ncol(inputWiggles) == 2 ) {
          	data <- inputWiggles
		names(data) <- c( "position", deparse(substitute(inputWiggles)) )
           } else { stop("Incorrect input format. Check documentation for details.
If this is a merged wiggle file, ensure first column is called 'position' and then run 'Compress'.") }
        }

	else if (is.vector(inputWiggles) && is.character(inputWiggles) ) {
	# allows unlimited number of files to compress on the same scale at the same time
	for (i in 1:length(inputWiggles)) {
	       	# to evaluate the strings extracted
		tmp = eval( parse( text = inputWiggles[i] ) )
		names(tmp) <- c( "position", inputWiggles[i] )
		if (i == 1) { tmp1 <- tmp }
		# use dplyr to join all data into a single data.frame for analysis
		else if (i == 2) {
		     data <- dplyr::full_join( tmp1, tmp, by="position" )
		} else {
		     data <- dplyr::full_join( data, tmp, by="position" )
		}
	}} else { stop("Incorrect input format. Check documentation for details.") }

	data <- data[ order( data$position ), ]	
	data2 <- Compress( data, window )
	
	message(paste0("Completed in ", round((proc.time()[3] - ptm[3])/60, 
        2), " min."))

	return(data2)
}

Compress <- function ( Data, window ) {

	if (!requireNamespace("dplyr", quietly = TRUE)) {
           stop("R package 'dplyr' needed for this function to work. Please install it.\n", 
            "install.packages('dplyr')", call. = FALSE)
    	}

       # determine size of output vector
       Size <- floor( ( max( Data$position ) - min( Data$position ) ) / window )
       
       out <- as.data.frame( matrix( nrow= Size, ncol= ncol(Data) ) )
       names(out) <- names(Data)
       
       # for progress bar
       pb <- txtProgressBar( min=0, max=Size, style=3 )
       
       # finds range of chr positions for each i
       # then determines range of rows with chr positions having values within this range
       # if statement deals with issues on chromosome ends or when entire region is missing from wiggle data
       # extracts rows and averages each column using DPLYR
       for ( i in 1:Size ) {
       	   Min <- min( Data$position ) + ( window * (i - 1) )
	   Max <- min( Data$position ) + ( window * i ) - 1
	   a <- min( which( Data$position >= Min ) )
	   b <- max( which( Data$position <= Max ) )
	   if ( a < b ) {
	      out[i,] <- dplyr::summarise_each( Data[a:b,], dplyr::funs(mean(., na.rm=TRUE )) )
	   }

	   # for progress bar
	   #Sys.sleep(0.1)
	   setTxtProgressBar( pb, i )
	}

	close(pb)
	
	out2 <- out[ which( out$position != 0 ),]
	return (out2)
}
