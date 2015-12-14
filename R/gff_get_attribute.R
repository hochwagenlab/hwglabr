#' Parse gff field
#'
#' This function allows you to parse the complex 'attributes' column of a gff file
#' and get only feature names or IDs, for example.
#' Adapted from function in the \href{https://bioconductor.org/packages/release/data/experiment/html/davidTiling.html}{davidTiling experimental package}.
#' Also posted on the Bioconductor support system \href{https://support.bioconductor.org/p/24657/}{here}.
#' @param gffColumn The column of the Gff file to parse, typically 'attributes'. No default.
#' @param field A string representing the field with the value you want
#' (formatted as 'field=value;'). No default.
#' @param attrSep A string representing the field separator (';' in the standard
#' gff format). Defaults to ';'.
#' @return A vector containing the specified field values of the same size as the
#' number of features in the gff file.
#' @section Details:
#' Typical use of this function will include first using gff_read(),
#' also in this package, to load gff file and then parsing the attributes field.
#' See examples below.
#' @examples
#' gff_read(s288C_annotation_R64_modified.gff)
#' 
#' gff_get_attribute(gff$attributes, 'Name')
#' 
#' gff_get_attribute(gff$attributes, 'ID')
#' @export

gff_get_attribute <- function(gffColumn, field, attrSep = ";") {
  s = strsplit(gffColumn, split = attrSep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    } else {
      rv = as.character(NA)
    }
    return(rv)
    })
}