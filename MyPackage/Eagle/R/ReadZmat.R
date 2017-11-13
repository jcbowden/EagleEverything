#' @title Read Z matrix
#' @description Read in the Z matrix that assigns groups/strains/lines to their trait measurements.
#' @param filename contains the name of the Z matrix file. The file name needs to be in quotes.
#' @details
#' The underlying linear mixed model is of the form
#' \deqn{Y = X \beta  + Z u_g + e}
#' where Z is a (n x n_g) matrix that contains ones and zeros, n is the number of trait measurements, and n_g 
#' is the number of groups/strains/lines. If n and n_g are the same, then there is no need to specify Z.   
#
#' 
#' A space separated text file is assumed.  Each row of the matrix contains multiple zeroes but only a  single one. 
#' The file cannot contain column or row headings. 
#'
#' 
#'
#'
#' @seealso \code{\link{ReadMarker}} and \code{\link{ReadPheno}}.
#' @return 
#' a data matrix is returned of the Z matrix.
#'
#' @examples
#' # Read in  example Z matrix from ./extdata/
#' 
#' # find the full location of the Z matrix data 
#' complete.name <- system.file('extdata', 'Z.txt', package='Eagle')
#'   
#' # read in Z matrix data 
#' Z_obj <- ReadZmat(filename=complete.name) 
#'                                
#'# look at first few rows of the Z matrix file
#' head(Z_obj)
#'
#'
ReadZmat  <- function( filename = NULL)
{
 Zfile <- fullpath(filename)
 error.code <-  check.inputs(file_genotype=filename)
 if(error.code){
    message(" ReadZmat has terminated with errors.")
   return(NULL)
  }
  Z <- fread(Zfile)
  Z <- as.matrix(Z)

 # do some checks 
 # number of unique values
 if(length(unique(as.vector(Z)))!=2){
    message(" ")
    message(" ERROR: The Z matrix file contains values other than 0 and 1.")
    message("   ")
    message("        ReadZmat has terminated with errors.")
    message(" ")
    return(NULL)
 }

 # checking that every row just doesnt contain zero values
 if(any(rowSums(Z)==0)){
   indx <- which(rowSums(Z)==0)
   if(length(indx)>0){

      message("  ")
      message(" ERROR:  The rows ", indx, " in the Z matrix have only 0 values.")
      message("         Each row must contain a single 1 value. ")
      message("   ")
      message("        ReadZmat has terminated with errors.")
      message(" ")
      return(NULL)

   }
 }


 # a row must contain only a single 1 value
 if(any(rowSums(Z)!=1)){
   indx <- which(rowSums(Z)!=1)
   if (length(indx) > 0){
      message("  ")
      message(" ERROR:  The rows ", indx, " in the Z matrix are incorrect.")
      message("         A row can only contains 0s and a single 1. ")
      message("   ")
      message("        ReadZmat has terminated with errors.")
      message(" ")
      return(NULL)
    }
  }


message("\n\n Loading Z matrix file ... \n\n")
message("                    Summary of Z matrix File  \n")
message("                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
message(" File name:                   ",  Zfile, "\n")
message(" Number of rows:              ", nrow(Z), "\n")
message(" Number of columns:           ", ncol(Z), "\n")
message(" First 5 rows of the Z matrix file are \n")

if(nrow(Z) > 5){
  for(ii in 1:5){
  message(cat(paste(Z[ii,], sep=" ")))
  }
} else {
  for(ii in 1:nrow(Z) ){
  message(cat(paste(Z[ii,], sep=" ")))
  }
}


message("\n\n")

return(Z)

}

