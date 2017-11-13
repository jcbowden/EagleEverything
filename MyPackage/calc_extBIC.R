 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, Zmat=NULL, quiet=TRUE)
 {
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, Z=Zmat, llim=-100,ulim=100)
   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

   extBIC <- BIC + 2 * lchoose(geno$dim_of_ascii_M[2], ncol(currentX) - 1)

    return(extBIC)
 }


