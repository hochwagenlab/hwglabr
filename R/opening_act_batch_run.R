#' Batch run of \code{\link{opening_act}} function
#'
#' This function will take a CSV file containing the path to multiple wiggle data folders
#' and argument labels and run \code{\link{opening_act}} function for each data set in the
#' file. It allows you to run a batch of data sets autmatically.
#' @param inputDataFile String indicating path to file containing instructions.
#' The file must be in CSV format with the following columns (column names in the file must
#' be reproduced exactly as below):
#' \enumerate{
#'   \item \strong{path} Path to wiggle data.
#'   \item \strong{relevantGenotype} String indicating the relevant strain mutations, to
#'   be used as the argument to \code{\link{opening_act}} of the same name.
#'   \item \strong{chipTarget} String indicating the ChIP target protein, to be used as
#'   the argument to \code{\link{opening_act}} of the same name.
#'   \item \strong{sampleID} The sample ID, to be used as the argument to
#'   \code{\link{opening_act}} of the same name.
#' }
#' No default.
#' @param useReadr Boolean indicating whether to use the much faster \code{read_tsv()} from
#' Hadley Wickham's readr package instead of base R's \code{read.table} when reading the wiggle
#' data. Used as the argument to \code{\link{readall_tab}} of the same name.
#' Defaults to \code{TRUE}.
#' @param runMetaORF Boolean indicating whether to run the meta ORF analysis. This analysis
#' typically takes about 30 minutes to run, so it may be useful to exclude it. Used as the
#' argument to \code{\link{opening_act}} of the same name. Defaults to \code{TRUE}.
#' @return The output of \code{\link{opening_act}} for each of the data sets included in the
#' batch. Check the documentation of \code{\link{opening_act}} for more detail.
#' @param userInput Boolean indicating whether to ask user to check the format of the
#' provided labels. Defaults to \code{FALSE}.
#' @examples
#' \dontrun{
#' opening_act_batch_run(inputDataFile="~/Desktop/inputs.csv")
#' opening_act_batch_run(inputDataFile="~/Desktop/inputs.csv", useReadr = FALSE,
#'                       runMetaORF = TRUE, userInput = FALSE)
#' }
#' @export

opening_act_batch_run <- function(inputDataFile, useReadr = TRUE, runMetaORF = TRUE,
                                  userInput = FALSE){
  ptm <- proc.time()
  
  # Have users check the arguments file for mistakes:
  if(userInput){
    # Ask user to make sure they provided valid arguments for opening_act
    title <- paste0('The arguments in the batch file will be used to name the final output folders.
For example, the "sampleID" label should identify the yeast strain, date,
and read mapping conditions, as in:
"AH119C-040114-sacCer3-2mis".\n
Please make sure all your provided labels are correct.')
    choices = c('Need to make changes, let me go back.',
                'Looking good, continue analysis!')
    answer <- menu(choices, graphics = FALSE, title)
    
    if(answer == 0 | answer == 1){
      stop('You chose to stop the function.', call. = FALSE)
    }
  }
  
  # Read in file with path to every wiggle data set and all annotations
  args_file <- read.csv(inputDataFile, stringsAsFactors = F)

  # Are all "paths" accessible?
  if (any(!file.exists(args_file[, 'path']))) {
    stop('Cannot seem to find the files for sample(s): ',
         paste(which(!file.exists(args_file[, 'path'])), collapse=", "), '\n',
         'Please check the provided path(s) to files.', call. = FALSE)
  }
  
  # Are all "paths" unique?
  if (length(unique(args_file[, 'path'])) != nrow(args_file)) {
    stop('Not all provided paths to data are unique. Please check the provided paths.',
         call. = FALSE)
  }
  
  # For each data set, read in wiggle data and run 'opening_act'
  for(i in 1:nrow(args_file)){
    wiggleData <- hwglabr::readall_tab(args_file[i, 'path'], localCopy = TRUE,
                                       useReadr = useReadr)
    
    opening_act(wiggleData, relevantGenotype = args_file[i, 'relevantGenotype'],
                chipTarget = args_file[i, 'chipTarget'],
                sampleID = args_file[i, 'sampleID'],
                userInput = FALSE, runMetaORF = runMetaORF)
    
    message('---')
    message('---')
    message('Ran "opening_act" on all data found in ', args_file[i, 'path'])
  }
  
  message()
  message('----------------------------------------')
  message('Ran "opening_act" on all ', nrow(args_file), ' datasets.')
  
  elapsed_time <- round((proc.time()[3] - ptm[3]), 1)
  if(elapsed_time < 60){
    message('\n...\nCompleted in ', elapsed_time, ' sec.')
  } else if(elapsed_time >= 60 & elapsed_time < 3600){
    message('\n...\nCompleted in ', round(elapsed_time / 60, 1), ' min.') 
  } else message('\n...\nCompleted in ', round(elapsed_time / 60 / 60, 1), ' h.')
  message('----------------------------------------')
}