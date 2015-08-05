## summarize findings from multiple_locus_am
## produce p values and amount of variance explained. 


#' @title Summary of multiple locus association mapping results
#' @description A summary of the results from \code{\link{multiple_locus_am}}
#' @param resam  the (list) object obtained from running \code{\link{multiple_locus_am}}.
#' @param  pheno  the (data frame) object  obtained  from running \code{\link{read.phenotypes}}.
#' @param geno   the (list) object obtained from running \code{\link{read.genotypes}}.
#' @param map   the (data frame) object obtained from running \code{\link{read.map}}.  
#' If not specifed, a generic map will be assumed. 
#' @param workingmemGb a numeric value. It specifies the amount of memory (in Gigabytes) available for reading analysis. 
#' @param numcores a numeric value for the number of cores that are available for distributed computing. 
#' @param  verbose a logical value. When \code{TRUE}, extra output is returned to the screen.
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
summarymlam <- function(resam=NULL, pheno=NULL, geno=NULL, map=NULL, workingmemGb, numcores, verbose)
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
  currentX <- .build_design_matrix(pheno=pheno, geno=geno, indxNA=resam$indxNA, colname.feffects=colname.feffects,
                                   verbose=verbose)


  ## add genetic marker effects 
  for(loc in am_res$Indx){
        currentX <- constructX(currentX=currentX, loci_indx=loc,
                               bin_path = resam$bin_path, 
                               dim_of_bin_M=geno[["dim_of_bin_M"]],
                               indxNA = resam$indxNA, map=map)
   }  ## end for loc

 ## calculate MMt
 MMt <- .calcMMt(geno, resam$workingmemGb, resam$numcores, resam$selected_loci, resam$verbose, resam$indxNA)

 ## calculate variance components of LMM
 eR <- emma.REMLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100)

 ## calculating p values of fixed marker effecs via Wald statistic
 pval <- vector("numeric", length(resam$Mrk))
 for(ii in resam$Mrk){
    indx <- which(colnames(currentX)==ii)
    L <- matrix(data=rep(0, ncol(currentX)), byrow=TRUE, nrow=1)
    L[indx] <- 1
    H <-  eR$vg * MMt + eR$ve * diag(1, nrow(MMt))
    Hinv <- try(solve(H))
    beta <- try(solve( t(currentX) %*% Hinv %*% currentX) %*% t(currentX) %*% Hinv %*% matrix(data=trait ,ncol=1)   )
    W <- t(L %*% beta) %*%
            solve( L %*% solve(t(currentX) %*% Hinv %*% currentX) %*% t(L) ) %*%
            (L %*% beta)
   pval[which(ii==resam$Mrk)] <- 1 - pchisq(W, 1) 
 }  ## end for ii in resam$Mrk

 ## print Annova table of results

 cat("Table of Results:  Marker loci in significant association with the trait \n")
 
if(length(selected_loci)==1 $  any(is.na(selected_loci)))
{
   cat(" No significant marker-trait associations have been found. \n\n")
} else {
   cat(sprintf("%15s     %10s        %10s        %10s    &10s    %10s \n",
        "Mrk Name", "Chrm", "Map Pos", "Col Number", "Additive effect", "p-value"))
   for(ii in resam$selected_loci)
   {
      beta_indx <- which(colnames(currentX) == map[[1]][ii])
      pval_indx <- which(ii==resam$selected_loci)
      if(length(beta_indx))
         stop(" Internal error in summarymlam function. Marker name mismatch. Contact package author.")

      cat(sprintf("%15s     %10s        %10f        %9i   %10f    %10f\n",
         map[[1]][ii], map[[2]][ii], map[[3]][ii], ii, beta[beta_indx], pval[pval_indx]))
   }  ## end for ii
} ## if else


  

}
