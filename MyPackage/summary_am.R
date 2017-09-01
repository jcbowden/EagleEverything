
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
#' proportion of phenotypes variance explained by  the different multiple-locus models. Each row 
#' in this table is the proportion of phenotype variance after the marker locus has been added to the 
#' multiple locus model. Our calculations of variance explained are based on Sun et al. (2010).  
#' @references  Sun G., Zhu C., Kramer  MH., Yang S-S., et al. 2010. Variation explained in mixed model association 
#' mapping. Heredity 105, 330-340. 
#' @examples
#'  \dontrun{
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #---------------
#'   # read the map 
#'   #---------------
#'   #
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   #------------------
#'   # read marker data
#'   #------------------
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8 gigabytes
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   #----------------------
#'   # read phenotype data
#'   #-----------------------
#'
#'   # Read in a plain text file with data on a single trait and two fixed effects
#'   # The first row of the text file contains the column names y, cov1, and cov2. 
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   #-------------------------------------------------------
#'   # Perform multiple-locus genome-wide association mapping 
#'   #-------------------------------------------------------                   
#'   res <- AM(trait = 'y',
#'                            fformula=c("cov1 + cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
#'
#'   #-----------------------------------------
#'   # Produce additional summary information 
#'   #------------------------------------------
#'
#'   SummaryAM(AMobj=res, pheno=pheno_obj, geno=geno_obj, map=map_obj)
#'  }
#'
#' 
#' 
#' @seealso \code{\link{AM}}
#'
SummaryAM <- function(AMobj=NULL, pheno=NULL, geno=NULL, map=NULL)
{

 if(is.null(AMobj)){
    message(" SummaryAM function requires AMobj object to be specified.")
    return(NULL)
    }
 if(is.null(pheno)){
    message(" SummaryAM function requires pheno parameter to be specified.")
    return(NULL)
    }
 if(is.null(geno)){
    message(" SummaryAM function requires geno parameter to be specified.")
    return(NULL)
    }
 if(!is.list(AMobj)){
    message(" SummaryAM function requires AMobj object to be a list object.")
    return(NULL)
   }
 if(!is.data.frame(pheno)){
    message(" SummaryAM function requires pheno object to be a data.frame object.")
    return(NULL)
    }
 if(!is.list(geno)){
   message(" SummaryAM function requires geno object to be a list object.")
    return(NULL)
   }

 if(is.null(map)){
   if(!AMobj$quiet ){
     message(" Map file has not been supplied. An artificial map is being created but this map is not used in the analysis. \n")
     message(" It is only used for the reporting of results. \n")
   }
   ## map has not been supplied. Create own map
   map <- data.frame(SNP=paste("M", 1:geno[["dim_of_ascii_M"]][2], sep=""), 
                     Chr=rep(1, geno[["dim_of_ascii_M"]][2]), 
                     Pos=1:geno[["dim_of_ascii_M"]][2])
  }

 ## check to make sure that null model is not being supplied
 if (length(AMobj$Mrk)==1){
   message(" No significant marker-trait associations have been found by AM. \n")
   message(" Nothing to summarize. \n")
   return()
 }

  ## build environmental effects design matrix
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

 ## calculating p values of fixed marker effect via Wald statistic
 mrks <- AMobj$Mrk[-1]  ## its -1 to remove the NA for the null model 
 pval <- vector("numeric", length(colnames(fullX)) )

 H <-  eR$vg * MMt + eR$ve * diag(1, nrow(MMt))
 Hinv <- try(solve(H))
 beta <- try(solve( t(fullX) %*% Hinv %*% fullX) %*% t(fullX) %*% Hinv %*% matrix(data=AMobj$trait ,ncol=1)   )
 df_pvalue <- NULL
 for(ii in colnames(fullX)  ){
    indx <- which(colnames(fullX)==ii)
    L <- matrix(data=rep(0, ncol(fullX)), byrow=TRUE, nrow=1)
    L[indx] <- 1
    W <- t(L %*% beta) %*%
            solve( L %*% solve(t(fullX) %*% Hinv %*% fullX) %*% t(L) ) %*%
            (L %*% beta)
   pval[which(ii==colnames(fullX)) ] <- 1 - pchisq(W, 1) 
 }  ## end for ii in AMobj$Mrk
df_pvalue <- data.frame("effects"=colnames(fullX), "p-value"=pval)
 ## print Annova table of results


  message(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
message("     Size and Significance of Effects in Final Model    \n")
  message(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")

  message(sprintf("%15s  %10s  %10s \n", "Name", "Additive effect", "p-value"))
  for(ii in colnames(fullX) )
  {
      indx <- which(colnames(fullX) == ii)
      message(sprintf("%15s  %10f         %.3E\n",
         ii, beta[indx], pval[indx ]))
  }  ## end for ii
 message("\n\n\n")
df_size <- data.frame("effect_names"=colnames(fullX), "estimate" = beta, "p_value"=pval)



 ##----------------------------------------------------------------------- 
 ## Variance explained - based on Sun et al. (2010). Heredity 105:333-340
 ##----------------------------------------------------------------------- 

 MMt <- MMt/max(MMt) + 0.05 * diag(nrow(MMt))  
 # base model
 basemod <- emma.MLE(y=AMobj$trait, X=baseX, K=MMt, llim=-100,ulim=100)
 base_logML <- basemod$ML

 # full model
  df_R <- NULL
  fullX <- baseX
  message(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
  message(" Proportion of Phenotype Variance Explained by Multiple-locus \n")
  message("             Association Mapping Model \n")
  message("  Marker loci which were found by AM() are added one at a time    \n")
  message(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
  message(sprintf("   %15s      %10s \n", "Marker name", "Proportion"))
  for(loc in AMobj$Indx[-1]){
        fullX <- constructX(fnameM=geno[["asciifileM"]],
                                currentX=fullX, loci_indx=loc,
                               dim_of_ascii_M=geno[["dim_of_ascii_M"]],
                               map=map)
        fullmod <- emma.MLE(y=AMobj$trait, X=fullX, K=MMt, llim=-100,ulim=100)
        full_logML <- fullmod$ML
        Rsq <- 1 - exp(-2/nrow(MMt) * (full_logML - base_logML))
        message(sprintf("  %+15s          %.3f\n",  paste("+ ",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])), Rsq))
        df_R <- rbind.data.frame(df_R, data.frame("Marker_name"=paste("+",as.character(AMobj$Mrk[which(loc==AMobj$Indx)])),
                                                  "Prop_var_explained"=Rsq))
   }  ## end for loc


  res <- list()
  res[["pvalue"]] <- df_pvalue
  res[["size"]] <- df_size
  res[["R"]] <- df_R

  
  invisible(res)
}
