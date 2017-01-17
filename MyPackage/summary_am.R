
GenomicRel = function(M){
#       markers must be 0, 1, 2 for homozygous, heterozygous and other homozygous

    M= M+1  ## since M is -1,0,1
    p1=round((apply(M,2,sum)+nrow(M))/(nrow(M)*2),3)
    p=2*(p1-.5)
    P = matrix(p,byrow=T,nrow=nrow(M),ncol=ncol(M))
    Z = as.matrix(M-P)

    b=1-p1
    c=p1*b
    d=2*(sum(c))

    ZZt = Z %*% t(Z)
    G = (ZZt/d)
    #invG=solve(G)

    return(G)
}




#' @title Summary of multiple locus association mapping results
#' @description    A summary function that provides additional information on the significant 
#'     marker-trait associations found by \code{\link{AM}}
#' @param  AMobj  the (list) object obtained from running \code{\link{AM}}. Must be specified. 
#' @param  pheno  the (data frame) object  obtained  from running \code{\link{ReadPheno}}. Must be specified. 
#' @param geno   the (list) object obtained from running \code{\link{ReadMarker}}. Must be specified. 
#' @param map   the (data frame) object obtained from running \code{\link{ReadMap}}. The default is to assume 
#'              a map object has not been supplied.   Optional.

#' @details
#'
#' \code{SummaryAM} produces two tables of results. First, a table of results is produced with 
#' the additive effect size and p-value for each 
#' fixed effect in the final model.  Second, a table of results is produced with the 
#' proportion of phenotypic variance explained by  the different multiple-locus models. Each row 
#' in this table is the proportion of phenotypic variance after the marker locus has been added to the 
#' multiple locus model. Our calculations of variance explained are based on Sun et al. (2010).  
#' @references  Sun G., Zhu C., Kramer  MH., Yang S-S., et al. 2010. Variation explained in mixed model association 
#' mapping. Heredity 105, 330-340. 
#' @examples
#'   #---------------
#'   # read the map 
#'   #---------------
#'   #
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file("extdata", "map.txt", 
#'                                    package="AMplus")
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   #------------------
#'   # read marker data
#'   #------------------
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8Gbytes
#'   complete.name <- system.file("extdata", "geno.ped", 
#'                                      package="AMplus")
#'   geno_obj <- ReadMarker(filename=complete.name,  type="PLINK", availmemGb=8) 
#'  
#'   #----------------------
#'   # read phenotypic data
#'   #-----------------------
#'
#'   # Read in a plain text file with data on a single trait and two covariates
#'   # The first row of the text file contains the column names "trait", "cov1", and "cov2". 
#'   complete.name <- system.file("extdata", "pheno.txt", package="AMplus")
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   #-------------------------------------------------------
#'   # Perform multiple-locus genome-wide association mapping 
#'   #-------------------------------------------------------                   
#'   res <- AM(trait = "trait",
#'                            fformula = c("cov1 + cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
#'
#'   #-----------------------------------------
#'   # Produce additional summary information 
#'   #------------------------------------------
#'
#'   SummaryAM(AMobj=res, pheno=pheno_obj, geno=geno_obj, map=map_obj)
#'
#'
#' 
#' @seealso \code{\link{AM}}
#'
SummaryAM <- function(AMobj=NULL, pheno=NULL, geno=NULL, map=NULL)
{

 if(is.null(AMobj))
    stop(" SummaryAM function requires AMobj object to be specified.", call. = FALSE)
 if(is.null(pheno))
    stop(" SummaryAM function requires pheno parameter to be specified.", call. = FALSE)
 if(is.null(geno))
    stop(" SummaryAM function requires geno parameter to be specified.", call. = FALSE)

 if(!is.list(AMobj))
    stop(" SummaryAM function requires AMobj object to be a list object.", call. = FALSE)
 if(!is.data.frame(pheno))
    stop(" SummaryAM function requires pheno object to be a data.frame object.", call. = FALSE)
 if(!is.list(geno))
    stop(" SummaryAM function requires geno object to be a list object.", call. = FALSE)

 if(is.null(map)){
   if(AMobj$quiet > 0){
     cat(" Map file has not been supplied. An artifical map is being created but this map is not used in the analysis. \n")
     cat(" It is only used for the reporting of results. \n")
   }
   ## map has not been supplied. Create own map
   map <- data.frame(SNP=paste("M", 1:geno[["dim_of_ascii_M"]][2], sep=""), 
                     Chr=rep(1, geno[["dim_of_ascii_M"]][2]), 
                     Pos=1:geno[["dim_of_ascii_M"]][2])
  }

 ## check to make sure that null model is not being supplied
 if (length(AMobj$Mrk)==1){
   cat(" No significant marker-trait associations have been found by AM. \n")
   cat(" Nothing to summarize. \n")
   return()
 }

  ## build enviornmental effects design matrix
  baseX <- .build_design_matrix(pheno=pheno,  indxNA=AMobj$indxNA, 
                                    fformula=AMobj$fformula,
                                   quiet=AMobj$quiet)


  ## add genetic marker effects 
  fullX <- baseX
  for(loc in AMobj$Indx){
        fullX <- constructX(fnameM=geno[["asciifileM"]], 
                              currentX=fullX, loci_indx=loc,
                               dim_of_ascii_M=geno[["dim_of_ascii_M"]],
                                map=map)
   }  ## end for loc

 ## calculate MMt
 MMt <- .calcMMt(geno, AMobj$availmemGb, AMobj$ncpu, AMobj$Indx, AMobj$quiet)

 ## calculate variance components of LMM
 eR <- emma.REMLE(y=AMobj$trait, X= fullX , K=MMt, llim=-100,ulim=100)

 ## calculating p values of fixed marker effecs via Wald statistic
 mrks <- AMobj$Mrk[-1]  ## its -1 to remove the NA for the null model 
 pval <- vector("numeric", length(colnames(fullX)) )

 H <-  eR$vg * MMt + eR$ve * diag(1, nrow(MMt))
 Hinv <- try(solve(H))
 beta <- try(solve( t(fullX) %*% Hinv %*% fullX) %*% t(fullX) %*% Hinv %*% matrix(data=AMobj$trait ,ncol=1)   )
 for(ii in colnames(fullX)  ){
    indx <- which(colnames(fullX)==ii)
    L <- matrix(data=rep(0, ncol(fullX)), byrow=TRUE, nrow=1)
    L[indx] <- 1
    W <- t(L %*% beta) %*%
            solve( L %*% solve(t(fullX) %*% Hinv %*% fullX) %*% t(L) ) %*%
            (L %*% beta)
   pval[which(ii==colnames(fullX)) ] <- 1 - pchisq(W, 1) 
 }  ## end for ii in AMobj$Mrk

 ## print Annova table of results


  cat(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
cat("     Size and Significance of Effects in Final Model    \n")
  cat(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

  cat(sprintf("%15s  %10s  %10s \n", "Name", "Additive effect", "p-value"))
  for(ii in colnames(fullX) )
  {
      indx <- which(colnames(fullX) == ii)
      cat(sprintf("%15s  %10f         %.3E\n",
         ii, beta[indx], pval[indx ]))
  }  ## end for ii
 cat("\n\n\n")




 ##----------------------------------------------------------------------- 
 ## Variance explained - based on Sun et al. (2010). Heredity 105:333-340
 ##----------------------------------------------------------------------- 

 MMt <- MMt/max(MMt) + 0.05 * diag(nrow(MMt))  
 # base model
 basemod <- emma.MLE(y=AMobj$trait, X=baseX, K=MMt, llim=-100,ulim=100)
 base_logML <- basemod$ML

 # full model
  fullX <- baseX
  cat(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
  cat(" Proportion of Phenotypic Variance Explained by Multiple-locus \n")
  cat("             Association Mapping Model \n")
  cat("  Marker loci which were found by AM() are added one at a time    \n")
  cat(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
  cat(sprintf("   %15s      %10s \n", "Marker name", "Proportion"))
  for(loc in AMobj$Indx[-1]){
        fullX <- constructX(fnameM=geno[["asciifileM"]],
                                currentX=fullX, loci_indx=loc,
                               dim_of_ascii_M=geno[["dim_of_ascii_M"]],
                               map=map)
        fullmod <- emma.MLE(y=AMobj$trait, X=fullX, K=MMt, llim=-100,ulim=100)
        full_logML <- fullmod$ML
        Rsq <- 1 - exp(-2/nrow(MMt) * (full_logML - base_logML))
        cat(sprintf("  %+15s          %.3f\n",  paste("+",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])), Rsq))
   }  ## end for loc



  

}
