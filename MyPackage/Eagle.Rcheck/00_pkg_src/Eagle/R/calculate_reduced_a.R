calculate_reduced_a <- function(varG=NULL, P=NULL, MMtsqrt=NULL, y=NULL, quiet=TRUE, message=message)
{

  ## internal function to AM
  if( !(nrow(P) ==  length(y))){
    message(" Error:  there is a problem with the  dimensions of  P, and/or the vector y.")
    message("         They should  be of the dimension (n x n), and a vector of length n.")
    message(" The dimensions are: \n")
    message(" dim(P)      = ", dim(P), "\n")
    message(" length(y)   = ", length(y), "\n")
    return(NULL)

  }

 if(is.null(varG)){
   message(" VarG must be specified.")
   return(NULL)
   }

  if(is.null(P)){
   message(" P must be specified")
   return(NULL)
   }


  if(is.null(y)){
   message(" y must be specified")
   return(NULL)
   }

    a <- varG * MMtsqrt %*% P %*% y

return(a)

}




