#' Signal around all bedfile midpoints
#'
#' This function allows you to pull out the ChIP signal centered around summits or 
#' midpoints of bedfiles if start and stop values are more than one base apart. It
#' takes as input wiggle data as a list of 16 chromosomes (output of \code{readall_tab()})
#' and a bedfile determining the positions to extract. Use with \code{signal_average()}
#' to calculate the mean around all midpoints.
#' Function written by Tovah Markowitz.
#' @param inputData A list of 16 chr wiggle data (output of readall_tab). No default.
#' @param bedData A data frame of bed data to extract signal around. No default.
#' @param extension Number indicating the number of bases around the summit or midpoints.
#' Value will be added both upstream and downstream. Default is 1000 bp.
#' @param onlyComplete Many wiggle files are missing data at random positions across the
#' genome. Boolean indicates if you want to only keep regions where every base is 
#' included within the wiggle file. Default is \code{TRUE}.
#' @return An R (dplyr) data frame with four columms: chromosome name, position (relative 
#' to midpoint), signal, and bed line number within its individual chromosome.
#' @examples
#' summit_plot(WT,red1_summit_bed)
#'
#' summit_plot(WT, S288Ccen, extension=2e4, onlyComplete = FALSE)
#' @export

# library(hwglabr)
# read in wiggle file using readall_tab

#bedData <- read.table("BedFiles/AH6408I-042415-SacCer3-2mis-M2-Tagnorm_P15_summits.bed",stringsAsFactors=F)
#bedData_factor <- read.table("BedFiles/AH6408I-042415-SacCer3-2mis-M2-Tagnorm_P15_summits.bed")

#inputData <- Red1WT_A

# want to include ability to plot around midpoints of bed_file as well before giving to Luis
# also gff?

signal_at_summit <- function(inputData, bedData, extension = 1000, onlyComplete = TRUE) {
ptm <- proc.time()

if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("R package 'dplyr' needed for this function to work. Please install it.\n", 
            "install.packages('dplyr')", call. = FALSE)
    }
    
# Check kind of numerals used in your input data (list of chromosome wiggle data)
check_S288C <- any(grep('chrI.', names(inputData), fixed = TRUE))
check_SK1 <- any(grep('chr01.', names(inputData), fixed = TRUE))

# Create vectors of chromosome numbers
chrom_S288C <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII",
                 "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII",
                 "chrXIV", "chrXV", "chrXVI")

chrom_SK1 <- c('chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07',
               'chr08', 'chr09', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
               'chr15', 'chr16')

# ensure that bedData does not contain factors for the chromosome names
if ( length( levels(bedData[,1]) ) != 0 ) {
	stop ("Chromosome names must be characters not factors.")
}

# ensure bedData is a bed file
if (grepl("chr", bedData[1,1]) & is.numeric(bedData[1,2]) & is.numeric(bedData[1,3])) {
	# ensure bedData is a summits file, or turn it into one
	if (sum( (bedData[,3]- bedData[,2]) != 1) != 0) {
		for (i in 1:nrow(bedData)) {
			bedData[i,2] <- floor(mean(c(bedData[i,2],bedData[i,3])))
			bedData[i,3]
		}
	}
} else stop ("bedData is in the incorrect format for this function.")

# ensure bedData uses the same reference genome as wiggle file
if (check_S288C) {
	cat ("Wiggle data is S288C.\n")
	if ( !( any(bedData[,1] == "chrI") ) ) {
		stop ("bedData and wiggle_data do not use the same reference genome.")
	} else {  chrom_nums <- chrom_S288C }
} else if (check_SK1) {
	cat ("Wiggle data is SK1.\n")
	if ( !( any(bedData[,1] == "chr01") ) ) {
		stop ("bedData and wiggle_data do not use the same reference genome.")
	} else { chrom_nums <- chrom_SK1 }
} else stop('Did not recognize reference genome.')

# ensure bedData is a bed file
if (grepl("chr",bedData[1,1]) & is.numeric(bedData[1,2]) & is.numeric(bedData[1,3])) {
	# ensure bedData is a summits file, or turn it into one
	if (sum( (bedData[,3]- bedData[,2]) != 1) != 0) {2
		# determine regions of genome to extract wiggle plot information
		size <- nrow (bedData)
		summit_range <- data.frame (chr=character(size), summit=numeric(size), start=numeric(size), 
			end=numeric(size), stringsAsFactors=F)
		for (i in 1:size) {
			summit_range$chr[i] 		<- bedData[i,1]
			summit_range$summit[i] 	<- floor(mean(c(bedData[i,2],bedData[i,3])))
			summit_range$start[i] 	<- summit_range$summit[i] - extension
			summit_range$end[i] 		<- summit_range$summit[i] + extension
		}
	} else {
		# determine regions of genome to extract wiggle plot information
		size <- nrow (bedData)
		summit_range <- data.frame (chr=character(size), summit=numeric(size), start=numeric(size), 
			end=numeric(size), stringsAsFactors=F)
		for (i in 1:size) {
			summit_range$chr[i]		<- bedData[i,1]
			summit_range$summit[i] 	<- bedData[i,2]
			summit_range$start[i] 	<- summit_range$summit[i] - extension
			summit_range$end[i]		<- summit_range$summit[i] + extension
		}
	}
} else stop ("bedData is in the incorrect format for this function.")

# determine regions of genome to extract wiggle plot information
size <- nrow (bedData)
summit_range <- data.frame (chr=character(size), summit=numeric(size), start=numeric(size), 
	end=numeric(size), stringsAsFactors=F)
for (i in 1:size) {
	summit_range$chr[i] 	<- bedData[i,1]
	summit_range$summit[i] 	<- bedData[i,2]
	summit_range$start[i] 	<- summit_range$summit[i] - extension
	summit_range$end[i] 	<- summit_range$summit[i] + extension
}

summit_data <- data.frame()
# iterate through list of chromosomes
for (i in 1:length(chrom_nums)) {
	# get wiggle information and rows from bedData associated with particular chromosome
	list_index <- grep(paste0(chrom_nums[i], '.'), names(inputData), fixed = TRUE)
	wiggle_data_chrom <- inputData[[list_index]]
	summit_range_chrom <- summit_range[summit_range[, 'chr'] == chrom_nums[i], ]
	cat (paste0(chrom_nums[i], ":\t", nrow(summit_range_chrom), " summits,\t"))
	count = 0
	for (j in 1:nrow(summit_range_chrom)) {
		# find start position or closest thing
		if ( length( which( wiggle_data_chrom[,1] >= summit_range_chrom[j,3]) )  == 0) {
			a <- 1
		} else {
			a <- min( which( wiggle_data_chrom[,1] >= summit_range_chrom[j,3]) )
		}
		
		# find end position or closest thing
		if ( length( which( wiggle_data_chrom[,1] <= summit_range_chrom[j,4]) ) == 0 ) {
			b <- nrow(wiggle_data_chrom)
		} else {
			b <- max( which( wiggle_data_chrom[,1] <= summit_range_chrom[j,4]) )
		}
		if (onlyComplete) {
			if ( (b-a) == (2 * extension) ) {
				count <- count + 1
				tmp <- wiggle_data_chrom[a:b,]
				colnames(tmp) <- c("position","signal")
				tmp[,1] <- tmp[,1] - summit_range_chrom[j,2]
				tmp2 <- data.frame( chr=rep( chrom_nums[i], nrow(tmp) ), tmp, line=rep( j, nrow(tmp) ), 
					stringsAsFactors=F)
				summit_data <- dplyr::bind_rows(summit_data,tmp2)
			}
		} else {
			# make sure I didn't end up with the entire chromosome
			if ( (a==1) & (b==nrow(wiggle_data_chrom)) ) {
				next
			} else if ( (a < b) & ((b-a) <= (2 * extension) ) ) { 
				count <- count + 1
				tmp <- wiggle_data_chrom[a:b,]
				colnames(tmp) <- c("position","signal")
				tmp[,1] <- tmp[,1] - summit_range_chrom[j,2]
				tmp2 <- data.frame( chr=rep( chrom_nums[i], nrow(tmp) ), tmp, line=rep( j, nrow(tmp) ), 
					stringsAsFactors=F)
				summit_data <- dplyr::bind_rows(summit_data,tmp2)
			}
		}
	}
	cat(paste0(count," summits mapped\n"))	
}

cat(paste0("Completed in ", round((proc.time()[3] - ptm[3])/60, 
        2), " min.\n"))
        
return(summit_data)

}