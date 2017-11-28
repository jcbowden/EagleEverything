check.for.NA.in.trait <- function(trait=NULL)
{
     ## internal function for AM 
     ## to return the positions of NA in a trait
     ## ordered for largest to smallest (this is important for ReshapeM_rcpp code

       ## check for NA's in trait
        indxNA <- which(is.na(trait))
        if(length(indxNA)==0){
          indxNA <- vector("numeric", 0)
        } else {
          ## place in reverse order
          indxNA <- sort(indxNA, decreasing = TRUE)
message(cat("\n\n WARNING!!!! The individuals in rows ", indxNA, " either have missing trait data "))
message("             and/or missing explanatory variable values. These individuals have ")
message(cat("             been removed from the analysis.  \n"))
          if(any(is.na(indxNA))){
            message("Error:  (internal).  indxNA contains NA values. ")
            message(" AM has terminated with errors. ")
            return(NULL)
          }
        }

      return(indxNA)
   }


