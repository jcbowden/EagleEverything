constructX <- function(fnameM=NULL, currentX=NULL, loci_indx=NULL,
                       availmemGb=8, dim_of_ascii_M=NULL,
                        map=NULL)
  {
    ## internal function for AM
    ## R function to construct the design matrix X
    ## Args
    ##   currentX    current model matrix
    ##   loci        the marker loci to be included as fixed QTL effects (additive model)

   if(is.na(loci_indx))
   {
     return(currentX)
   } else {
       genodat <- extract_geno(fnameM=fnameM, colnum=loci_indx,
                           availmemGb=availmemGb, dim_of_ascii_M=dim_of_ascii_M)
      newX <- cbind(currentX, genodat)
      colnames(newX) <- c(colnames(currentX), as.character(map[[1]][loci_indx])) ## adding col names to new X  
      return(newX)
   }
  }



