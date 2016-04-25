
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


## summarize findings from multiple_locus_am
## produce p values and amount of variance explained. 


#' @title Summary of multiple locus association mapping results
#' @description A summary of the results from \code{\link{multiple_locus_am}}
#' @param resam  the (list) object obtained from running \code{\link{multiple_locus_am}}.
#' @param  pheno  the (data frame) object  obtained  from running \code{\link{read.phenotypes}}.
#' @param geno   the (list) object obtained from running \code{\link{read.genotypes}}.
#' @param map   the (data frame) object obtained from running \code{\link{read.map}}.  
#' If not specifed, a generic map will be assumed. 
#' @details
#' \code{summarymlam} produces a summary table of the results from running \code{\link{multiple_locus_am}}. The 
#' summary table contains a separate row for each marker locus that is significantly associated with the trait.
#' Each row of the table contains the marker name, which chromosome it resides, its position, the size of 
#' the additive effect of the putative qtl that is in linkage disequilibrium with the marker locus, 
#' and its p-value in the linear mixed model.  The p-value is calculated from the Wald statistic.
#'
#' The amount of phenotypic variance explained by the selected loci is also reported.    
#'
#' @seealso \code{\link{multiple_locus_am}}
#'
summarymlam <- function(resam=NULL, pheno=NULL, geno=NULL, map=NULL)
{

 if(is.null(resam))
    stop(" summarymlam function requires resam object to be specified.")
 if(is.null(pheno))
    stop(" summarymlam function requires pheno parameter to be specified.")
 if(is.null(geno))
    stop(" summarymlam function requires geno parameter to be specified.")
 if(is.null(map))
    stop(" summarymlam function requires map parameter to be specified.")

 if(!is.list(resam))
    stop(" summarymlam function requires resam object to be a list object.")
 if(!is.data.frame(pheno))
    stop(" summarymlam function requires pheno object to be a data.frame object.")
 if(!is.list(geno))
    stop(" summarymlam function requires geno object to be a list object.")
 if(!is.data.frame(map))
   stop(" summarymlam function requires map to be a data.frame object.")



  ## build enviornmental effects design matrix
  baseX <- .build_design_matrix(pheno=pheno,  indxNA=resam$indxNA, 
                                    colname.feffects=resam$colname.feffects,
                                   verbose=resam$verbose)


  ## add genetic marker effects 
  fullX <- baseX
  for(loc in resam$Indx){
        fullX <- constructX(currentX=fullX, loci_indx=loc,
                               bin_path = resam$bin_path, 
                               dim_of_bin_M=geno[["dim_of_bin_M"]],
                               indxNA = resam$indxNA, map=map)
   }  ## end for loc

 ## calculate MMt
 MMt <- .calcMMt(geno, resam$workingmemGb, resam$numcores, resam$Indx, resam$verbose, resam$indxNA)

 ## calculate variance components of LMM
 eR <- emma.REMLE(y=resam$trait, X= fullX , K=MMt, llim=-100,ulim=100)

 ## calculating p values of fixed marker effecs via Wald statistic
 pval <- vector("numeric", length(resam$Mrk))
 for(ii in resam$Mrk){
    indx <- which(colnames(fullX)==ii)
    L <- matrix(data=rep(0, ncol(fullX)), byrow=TRUE, nrow=1)
    L[indx] <- 1
    H <-  eR$vg * MMt + eR$ve * diag(1, nrow(MMt))
    Hinv <- try(solve(H))
    beta <- try(solve( t(fullX) %*% Hinv %*% fullX) %*% t(fullX) %*% Hinv %*% matrix(data=resam$trait ,ncol=1)   )
    W <- t(L %*% beta) %*%
            solve( L %*% solve(t(fullX) %*% Hinv %*% fullX) %*% t(L) ) %*%
            (L %*% beta)
   pval[which(ii==resam$Mrk)] <- 1 - pchisq(W, 1) 
 }  ## end for ii in resam$Mrk

 ## print Annova table of results


cat("\n\n                       Summary of Association Mapping Results         \n\n")

 
if(length(resam$Indx)==1 &  any(is.na(resam$Indx)))
{
   cat(" No significant marker-trait associations have been found. \n\n")
} else {
   cat(sprintf("%15s  %10s      %10s   %15s   %10s  %10s \n",
        "Mrk Name", "Chrm", "Map Pos", "Col Number", "Additive effect", "p-value"))
   for(ii in resam$Indx)
   {
      beta_indx <- which(colnames(fullX) == map[[1]][ii])
      pval_indx <- which(ii==resam$Indx)
      if(length(beta_indx)==0)
         stop(" Internal error in summarymlam function. Marker name mismatch. Contact package author.")

      cat(sprintf("%15s  %10s      %10f     %9i   %10f         %.3E\n",
         map[[1]][ii], map[[2]][ii], map[[3]][ii], ii, beta[beta_indx], pval[pval_indx]))
   }  ## end for ii
} ## if else


 
 ## Variance explained - based on Sun et al. (2010). Heredity 105:333-340

 MMt <- MMt/max(MMt) + 0.05 * diag(nrow(MMt))  
 # base model
 basemod <- emma.MLE(y=resam$trait, X=baseX, K=MMt, llim=-100,ulim=100)
 base_logML <- basemod$ML

 # full model
 fullmod <- emma.MLE(y=resam$trait, X=fullX, K=MMt, llim=-100,ulim=100)
 full_logML <- fullmod$ML
 
 Rsq <- 1 - exp(-2/nrow(MMt) * (full_logML - base_logML))

 cat("\n\n       Amount of phenotypic variance explained is ", Rsq, " \n\n")
  

}
