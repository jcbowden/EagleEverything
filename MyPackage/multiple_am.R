.form_results <- function(trait, selected_loci, map, colname.trait, colname.feffects, bin_path, indxNA,
                           numcores, workingmemGb, verbose, herit, extBIC, stopping_rule)
{
  if (length(selected_loci) > 1){
#     if(stopping_rule=="BIC"){
#     sigres <- list(trait=trait, 
#                    colname.trait = colname.trait, 
#                    colname.feffects = colname.feffects,
#                    bin_path = bin_path,
#                    indxNA = indxNA,
#                    Mrk=map[[1]][selected_loci[-length(selected_loci)]], 
#                    Chr=map[[2]][selected_loci[-length(selected_loci)]], 
#                    Pos=map[[3]][selected_loci[-length(selected_loci)]], 
#                    Indx=selected_loci[-length(selected_loci)],
#                    numcores=numcores,
#                    workingmemGb=workingmemGb,
#                    verbose=verbose,
#                    herit=herit[selected_loci[-length(selected_loci)]],
#                    extBIC=extBIC[selected_loci[-length(selected_loci)]])
#      } else {
      sigres <- list(trait=trait,
                    colname.trait = colname.trait, 
                    colname.feffects = colname.feffects,
                    bin_path = bin_path,
                    indxNA = indxNA,
                    Mrk=map[[1]][selected_loci], 
                    Chr=map[[2]][selected_loci], 
                    Pos=map[[3]][selected_loci], 
                    Indx=selected_loci,
                    numcores=numcores,
                    workingmemGb=workingmemGb,
                    verbose=verbose,
                    herit=herit, 
                    extBIC=extBIC)





#     } ## end inner if
  } else {
      sigres <- NA
  }
 return(sigres)
}

.print_title <- function(){
    ## internal fuction: use only in multiple_locus_am function
    ## title
    cat("\n\n\n\n")
    cat("            Multiple Locus Association Mapping via WGAM\n")
    cat("                       Version 1.0 \n\n")
}


.build_design_matrix <- function(pheno=NULL, geno=NULL, indxNA=NULL, colname.feffects=NULL, verbose=FALSE  )
{
   ## internal fuction: use only in multiple_locus_am function and summaryam function
   ## build design matrix given character vector colname.feffects of column names

   ## assign model matrix X
   if(is.null(colname.feffects))
   {  ## trait + intercept being fitted only
      if(length(indxNA) > 0){
         Xmat <- matrix(data=1, nrow=nrow(pheno[-indxNA,]), ncol=1)

      } else {
        Xmat <- matrix(data=1, nrow=nrow(pheno), ncol=1)
      }
      colnames(Xmat) <- "intercept"
   } else {
      ## trait + fixed effects being fitted. 
     if(length(indxNA)==0)
     {
        mf <- paste(colname.feffects, collapse=" + ")
        mf <- paste(" ~ ", mf, sep="")
        mf <- as.formula(mf)
        Xmat <- model.matrix(mf, data=pheno)
     }  else {
        # there is an issue with creating Xmat when it includes
        # factors that have some of their levels removed. 
        ph <- pheno[-indxNA,]
        for(ii in colname.feffects){
           if(is.factor(ph[,ii])){
              ph[,ii] <- as.factor(as.character(ph[,ii]))
           }
        }  ## for    
        mf <- paste(colname.feffects, collapse=" + ")
        mf <- paste(" ~ ", mf, sep="")
        mf <- as.formula(mf)
        Xmat <- model.matrix(mf, data=ph)
     } ## if else (length(indxNA)==0)
   } 

 if (verbose){
   cat("Dimension of design matrix, before addition of marker fixed effects is ", nrow(Xmat), "rows and ", ncol(Xmat), "columns.\n") 
 }

if(!is.matrix(Xmat))
   Xmat <- matrix(data=Xmat, ncol=1)

  return(Xmat)
}


.calcMMt <- function(geno, workingmemGb, numcores, selected_loci, verbose, indxNA)
  {
    ## internal function: used only in multilocus_loci_am and summaryam
    ## values passed by environments
    MMt <- calculateMMt(geno=geno[["binfileM"]], workingmemGb=workingmemGb, 
                           numcores=numcores, 
                           dim_of_bin_M = geno[["dim_of_bin_M"]], 
                           selected_loci=selected_loci, verbose = verbose) 
    gc()

    if(length(indxNA)> 0 )
        MMt <- MMt[-indxNA, -indxNA]

    ## Trick for dealing with singular MMt due to colinearity
    MMt <- MMt/max(MMt) + diag(0.05, nrow(MMt)) 

    return(MMt)
  }

  .calcVC <- function(trait, currentX, MMt)
  {
    ## perform likelihood ratio test for variance component Var_g
    res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100)
    return(list("vg"=res_full$vg, "ve"=res_full$ve))

  }

 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, verbose=FALSE)
 {
   ## internal function: use in multiple_loucs_am only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100)
   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

   extBIC <- BIC + 2 * lchoose(geno$dim_of_bin_M[2], ncol(currentX) - 1)  

   if(verbose){
         cat(" new BIC = ", BIC, "\n")
         cat(" New extBIC = ", extBIC, "\n")
    }
    return(extBIC)
 }

 .print_header <- function(){
   cat(" No marker-trait associations found .... \n\n\n")
   cat("\n\n                           FINAL MODEL  \n")
   cat(" --------------------------------------------------------------------  \n")
 }

.print_findings <- function(selected_loci, map, herit, extBIC, stopping_rule)
{
    print(extBIC)

 #    if(stopping_rule=="BIC"){
 #        if (length(selected_loci[-length(selected_loci)]) == 1 & any(is.na(selected_loci)))
 #        {
 #            ## its -length() because while(continue) doesn't stop until it has gone 1 beyond
 #            ## the best set of selected_loci. 
 #             cat(" No significant marker-trait associations have been found. \n\n")
 #        }  else {
 #             cat(sprintf("%15s     %10s    %10s     %10s     %10s   %10s\n", 
 #                       "Mrk Name", "Chrm", "Map Pos", "Col Number", "Heritability", "extBIC"))
 #             ## its [-length()] because extBIC doesnt caouse the loop to stop until it has 
 #                 ## gone one beyond the best model.
 #              for(ii in selected_loci[-length(selected_loci)])
 #                     cat(sprintf("%15s     %10s        %10f        %9i    %4f   %4f \n", 
 #                         map[[1]][ii], map[[2]][ii], map[[3]][ii], ii, 
 #         herit[which(ii==selected_loci)], extBIC[which(ii==selected_loci)] ))
 #           cat(" --------------------------------------------------------------------  \n")
 #           cat("\n\n")
 #        }
 #  
 #   }  else {
       ### stopping_rule == "h"
              if (length(selected_loci) == 1 & any(is.na(selected_loci)))
              {
                   cat(" No significant marker-trait associations have been found. \n\n")
              }  else {
                   cat(sprintf("%15s     %10s        %10s        %10s       %10s   %10s\n",
                        "Mrk Name", "Chrm", "Map Pos", "Col Number", "Heritability", "extBIC"))

                   for(ii in selected_loci)
                           cat(sprintf("%15s     %10s        %10f        %9i    %4f  %4f\n",
                               map[[1]][ii], map[[2]][ii], map[[3]][ii], ii, 
                               herit[which(ii==selected_loci)], 
                               extBIC[which(ii==selected_loci)] ))
                 cat(" --------------------------------------------------------------------  \n")
                 cat("\n\n")
              }
 #  }


}  ## end function print.finals

 .print_results <- function(itnum, selected_loci, map,h, extBIC)
 { 
    cat(" Significant marker-trait association found ... \n")
    cat(" Results after iteration ", itnum, "\n")
    cat(sprintf("%15s     %10s        %10s        %10s      %10s  %10s\n", 
         "Mrk Name", "Chrm", "Map Pos", "Col Number","Heritability", "extBIC"))
    for(ii in selected_loci)
       cat(sprintf("%15s     %10s        %10f        %9i      %3f  %4f \n", 
        map[[1]][ii], map[[2]][ii], map[[3]][ii], ii, h[which(ii==selected_loci)],
        extBIC[which(ii==selected_loci)] ))
    cat("\n")
 }



  .find_qtl <- function(geno, workingmemGb, indxNA, selected_loci, MMt, invMMt, best_ve, best_vg, 
                       currentX, error_checking, numcores, verbose, trait)
  {
    ##  internal function: use only with multiple_locus_am
    if(verbose) cat(" Calculating H matrix. \n")
    H <- calculateH(MMt=MMt, varE=best_ve, varG=best_vg)
    if(verbose) cat(" Calculating P matrix. \n")
    P <- calculateP(H=H, X=currentX)
    if (verbose)
              cat(" Calculating  square root of M %*% t(M) and it's inverse. \n")
    MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking, numcores=numcores )

    if (verbose)
            cat(" Calculating BLUPs for dimension reduced model. \n")

    hat_a <- calculate_reduced_a(varG=best_vg, P=P, 
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                       y=trait, verbose = verbose )  

    if (verbose) 
             cat(" Calculating variance of BLUPs for dimension reduced model. \n")
    var_hat_a    <- calculate_reduced_vara(X=currentX, varE=best_ve, varG=best_vg, invMMt=invMMt, 
                                                MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                                                verbose = verbose)
   
    if (verbose)
         cat(" Calculating BLUPs and their variances for full model. \n")
    bin_path <- dirname(geno[["binfileM"]])
    a_and_vara  <- calculate_a_and_vara(bin_path=bin_path,  maxmemGb=workingmemGb, 
                                            dims=geno[["dim_of_bin_M"]],
                                            selectedloci = selected_loci,
                                            invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]],
                                            transformed_a=hat_a, 
                                            transformed_vara=var_hat_a,
                                            indxNA = indxNA,
                                            verbose=verbose)
  
    ## outlier test statistic
    if (verbose) cat(" Calculating outlier test statistics. \n")
    tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
    indx <- which(tsq == max(tsq, na.rm=TRUE))   ## index of largest test statistic. However, need to account for other loci 
                                         ## already having been removed from M which affects the indexing

    ## taking first found qtl
    indx <- indx[1]

    orig_indx <- seq(1, geno[["dim_of_bin_M"]][2])  ## 1:ncols

    return(orig_indx[indx])
}

#' @title multiple locus association mapping 
#' @description \code{multiple_locus_am} is used to perform multiple locus 
#' association mapping via multi-locus whole-genome  association mapping (MWAM)
#' @param numcores a numeric value for the number of cores that are available for distributed computing. 
#' @param workingmemGb a numeric value. It specifies the amount of memory (in Gigabytes) available for reading analysis. 
#' @param colname.trait  the name of the column in the phenotypic data file that contains the trait data. The name is case sensitive. 
#' @param colname.feffects   a character vector containing the column names of 
#'                        the explanatory/independent variables in the phenotypic data file. If
#'                        not specified, only an overall mean will be fitted.
#' @param map   the (data frame) object obtained from running \code{\link{read.map}}. If not specifed, a generic map will 
#'              be assumed. 
#' @param stopping_rule  takes a character value of "BIC" or "h" (for heritability). See details.
#' @param pheno  the (data frame) object  obtained  from running \code{\link{read.phenotypes}}.
#' @param geno   the (list) object obtained from running \code{\link{read.genotypes}}.
#' @param alpha  the type 1 error rate where setting \code{alpha} to 0.05 say is a 5\% error rate.
#' @param  error_checking a logical value. When \code{TRUE}, 
#' the numericial stability of the dimension reduction is checked. That is, individuals 
#' with near identical marker genotypes can cause numerical issues, and are reported 
#' when \code{error_checking} has been set to \code{TRUE}. 
#' @param  verbose      a logical value. When \code{TRUE}, extra output is returned 
#'  to the screen for monitoring progress. 
#' @param maxit     an integer value for the maximum number of forward steps to be performed. That is, it is the maximum number of 
#' qtl that are to be included in the model. 
#' @details
#' The steps to running \code{multiple_locus_am} are as follows:
#' \itemize{
#' \item{Step 1:}{ Read in  genotypic information. Use \code{\link{read.genotypes}} to read in the marker data.}
#' \item{Step 2:}{ Read in phenotypic information.  Use \code{\link{read.phenotypes}} to read in the phenotypic data.}
#' \item{Step 3:}{ Read in map information. Use \code{\link{read.map}} to read in 
#'  the marker map. Omit this step if marker map is unknown. A generic map will be created.}
#' \item{Step 4:}{ Perform WGAM analysis. Use \code{multiple_locus_am} to identify multiple 
#' significant marker-trait associations, simultaneously.}  
#' }
#'
#' The genotypic file must not contain missing genotypes.
#'
#' The phenotypic file can contain missing information, coded as \code{NA}, but only 
#' in the trait data. The fixed effects data (i.e. the covariates and/or factors)
#' cannot contain any missing data.  The phenotypic file is allowed to contain multiple traits
#' and fixed effects.  
#'  
#' The trait data and fixed effects data are  specified by setting \code{colname.trait} 
#'  and \code{colname.feffects}, respectively. \code{colname.trait} can only contain a single
#'  character string for the trait name.  \code{colname.feffects} can be a character vector 
#'  if multiple fixed effects are to be included in the linear mixed model. Whether the 
#'  fixed effects are treated as a covariate or factor is dependent upon the class of 
#'  the associated data columns in the data frame obtained from \code{\link{read.phenotypes}}. 
#'
#'    STILL BEING WRITTEN .... 
#' The \code{multiple_locus_am} function is an R/Rcpp implementation of multi-locus whole-genome 
#'  association mapping. The method is a 
#' multiple locus method that is a hybrid of whole-genome and multi-locus association mapping. 
#' Briefly, a multiple locus model
#' is built iteatively, by fitting a whole-genome model at each step. It differs from 
#' whole-genome mapping methods because we
#' get a parmonious set of marker loci that are in strongest association with a trait.  Also, it differs from multi-locus 
#' association mapping methods because at each iteration of the model building process, we fit all snp loci 
#' simultaneously, as opposed to fitting them one at a time. 
#'
#' NEEDS MORE  Add something about stopping rule
#'
#'
#' @seealso \code{\link{read.genotypes}}, and \code{\link{read.phenotypes}}.
#'
#' @return
#' something here .... 
#' @examples
#'   # READ MAP INFORMATION
#'   map.file.loc <- system.file("extdata", "mapexample.txt", 
#'                                    package="MWAM")
#'   map.df <- read.map(path=dirname(map.file.loc),  
#'                       file_map=basename(map.file.loc)) 
#'
#'
#'   # READ GENOTYPE INFORMATION
#'   #  0,1 genotypes
#'   #  column wise marker data
#'   gen.file.loc <- system.file("extdata", "genoexampleCwise.txt", 
#'                                      package="MWAM")
#'   geno.list <- read.genotypes(path=dirname(gen.file.loc), 
#'                               columnwise=TRUE, AA=0, BB=1, 
#'                               file_genotype=basename(gen.file.loc),  
#'                               workingmemGb=4) 
#'  
#'   # READ PHENOTYPIC INFORMATION
#' phen.file.loc <- system.file("extdata", "phenoexample.csv", package="MWAM")
#'   
#' phenodf <- read.phenotypes(path=dirname(phen.file.loc),  
#'                              file_phenotype=basename(phen.file.loc), 
#'                              csv=TRUE) 
#'                                
#'  # PERFORM Whole-Genome Multi-locus Association Mapping
#'   res <- multiple_locus_am(colname.trait = "trait1",
#'                            colname.feffects = c("cov1","cov2", "fac"),
#'                            map = map.df,
#'                            pheno = phenodf,
#'                            geno = geno.list, verbose=TRUE, workingmemGb=4)
#'
#'
#'
#'
#'
multiple_locus_am <- function(numcores=1,workingmemGb=8, 
                              stopping_rule = c("h","BIC"),
                              colname.trait = NULL, 
                              colname.feffects  = NULL,
                              map = NULL,
                              pheno=NULL, 
                              geno=NULL, 
                              alpha=0.05, error_checking=FALSE, 
                              verbose=FALSE,
                              maxit=20){
 ## Core function for performing whole genome association mapping with EMMA
 ## Args
 ## numcores        number of cores available for computation
 ## memoryGb        maximum amount of working memory available for computation
 ## pheno           data frame 
 ##                 remaining columns are explanatory variables to include in the model. If a numeric vector, then it 
 ##                 is only a response to be fitted. 
 ## geno            if geno is a matrix or data frame, then the user has not read.genotypes and a bin packed file
 ##                 has not been created. If it is a character string, then it is the file location of the binary packed files. 
 ## alpha           significance level at which to perform the likelihood ratio test
 ## error_checking  when true, it performs some checks of the calculations
 ## maxit           maximum number of qtl to include in the model


 ## check parameter inputs
 stopping_rule <- match.arg(stopping_rule)
 check.inputs.mlam(numcores, workingmemGb, colname.trait, colname.feffects, 
                     map, pheno, geno, alpha )
 ## ADD CHECK TO MAKE SURE PHENO AND GENO NROWS ARE CORRECT - NOT BEING CHECKED AT THE MOMENT
 continue <- TRUE
 selected_loci <- NA
 herit <- NA
 extBIC <- NA
 ## assign trait 
 trait <-  pheno[[colname.trait]]

 ## check for NA's in trait
 indxNA <- check.for.NA.in.trait(trait=trait)

 ## remove missing observations from trait
 if(length(indxNA)>0){
    trait <- trait[-indxNA]
 }

 ## build design matrix currentX

Args <- list(pheno=pheno, geno=geno, indxNA=indxNA, colname.feffects=colname.feffects, verbose=verbose )
currentX <- do.call(.build_design_matrix, Args)


 ## print tile
 .print_title()

 ## Initialization
 old_extBIC <- 1e10 ; itnum <- 1
 h <- 1

 ## main loop

 while(continue){
       cat(" Performing iteration ... ", itnum, "\n")

    
    ## Calculate variance components
    ##---------------------------
    #  if(itnum == 1){
         ## Calculate MMt and its inverse
         Args <- list(geno=geno,workingmemGb=workingmemGb,
                    numcores=numcores,selected_loci=selected_loci,
                    verbose=verbose, indxNA=indxNA)
         MMt <- do.call(.calcMMt, Args)
         invMMt <- chol2inv(chol(MMt))
         gc()
         vc <- .calcVC(trait=trait, currentX=currentX,MMt=MMt)
         best_ve <- vc[["ve"]]
         best_vg <- vc[["vg"]]
    #  }
 
    ## Calculate extBIC
    ##---------------------
    new_extBIC <- .calc_extBIC(trait, currentX,MMt, geno, verbose)
    h <- best_vg/(best_vg + best_ve)

    cat(" old BIC ", old_extBIC, " new BIC ", new_extBIC, "h", h,  "\n")

    ## find QTL
    ##-----------
    old_extBIC <- new_extBIC
    ARgs <- list(geno=geno,workingmemGb=workingmemGb, indxNA=indxNA, selected_loci=selected_loci,
                     MMt=MMt, invMMt=invMMt, best_ve=best_ve, best_vg=best_vg, currentX=currentX,
                     error_checking=error_checking,
                      numcores=numcores, verbose=verbose, trait=trait )
    new_selected_locus <- do.call(.find_qtl, ARgs)  ## memory blowing up here !!!! 

    if(any(is.na(selected_loci))  & length(selected_loci)==1){
        selected_loci <- new_selected_locus
        herit <- h
        extBIC <- new_extBIC
    } else {
        selected_loci <- c(selected_loci, new_selected_locus)
        herit <- c(herit, h)
        extBIC <- c(extBIC, new_extBIC)
    }

    .print_results(itnum, selected_loci, map, herit, extBIC)

    currentX <- constructX(currentX=currentX, loci_indx=new_selected_locus,
                               bin_path = dirname(geno[["binfileM"]]),
                               dim_of_bin_M=geno[["dim_of_bin_M"]],
                               indxNA = indxNA,
                               map=map, workingmemGb = workingmemGb)

    ## check if loop needs to be stopped  
    if(stopping_rule == "BIC" && new_extBIC > old_extBIC){
        continue <- FALSE
        .print_header()
        .print_findings(selected_loci, map, herit, extBIC, stopping_rule)
    } else if (stopping_rule == "h" && h < 0.01){
       continue <- FALSE
       .print_header()
       .print_findings(selected_loci, map, herit, extBIC, stopping_rule)
    } else {
       # QTL present
        old_extBIC <- new_extBIC
    }  ## if new_extBIC 

     itnum <- itnum + 1
     if(itnum > maxit)
         continue <- FALSE 
 
  }  ## end while continue


  sigres <- .form_results(trait, selected_loci, map, colname.trait, colname.feffects, dirname(geno[["binfileM"]]), indxNA, numcores, workingmemGb, verbose, herit, extBIC, stopping_rule)   



return( sigres )

} ## end multiple_locus_am















