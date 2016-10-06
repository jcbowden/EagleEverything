.form_results <- function(trait, selected_loci, map,  feffects, indxNA,
                           ncpu, availmemGb, verbose, herit, extBIC )
{
  if (length(selected_loci) > 1){
   sigres <- list(trait=trait,
                    feffects = feffects,
                    indxNA = indxNA,
                    Mrk=map[[1]][selected_loci], 
                    Chr=map[[2]][selected_loci], 
                    Pos=map[[3]][selected_loci], 
                    Indx=selected_loci,
                    ncpu=ncpu,
                    availmemGb=availmemGb,
                    verbose=verbose,
                    herit=herit, 
                    extBIC=extBIC)
  } else {
   sigres <- list(trait=trait,
                    feffects = feffects,
                    indxNA = indxNA,
                    Mrk=NA,
                    Chr=NA,
                    Pos=NA,
                    Indx=selected_loci,
                    ncpu=ncpu,
                    availmemGb=availmemGb,
                    verbose=verbose,
                    herit=herit, 
                    extBIC=extBIC)
  }
return(sigres)
}

.print_title <- function(){
    ## internal fuction: use only in AM function
    ## title
    cat("\n\n\n\n")
    cat("            Multiple Locus Association Mapping via WGAM\n")
    cat("                       Version 1.0 \n\n")
}


.build_design_matrix <- function(pheno=NULL, geno=NULL, indxNA=NULL, feffects=NULL, verbose=FALSE  )
{
   ## internal fuction: use only in AM function and summaryam function
   ## build design matrix given character vector feffects of column names

   ## assign model matrix X
   if(is.null(feffects))
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
        mf <- paste(feffects, collapse=" + ")
        mf <- paste(" ~ ", mf, sep="")
        mf <- as.formula(mf)
        Xmat <- model.matrix(mf, data=pheno)
     }  else {
        # there is an issue with creating Xmat when it includes
        # factors that have some of their levels removed. 
        ph <- pheno[-indxNA,]
        for(ii in feffects){
           if(is.factor(ph[,ii])){
              ph[,ii] <- as.factor(as.character(ph[,ii]))
           }
        }  ## for    
        mf <- paste(feffects, collapse=" + ")
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


.calcMMt <- function(geno, availmemGb, ncpu, selected_loci, verbose, indxNA)
  {
    ## internal function: used only in multilocus_loci_am and summaryam
    ## values passed by environments
    cat(" Inside calcMMt ... \n")
    MMt <- calculateMMt(geno=geno[["binfileM"]], availmemGb=availmemGb, 
                           ncpu=ncpu, 
                           dim_of_bin_M = geno[["dim_of_bin_M"]], 
                           selected_loci=selected_loci, verbose = verbose) 
    gc()

    if(length(indxNA)> 0 )
        MMt <- MMt[-indxNA, -indxNA]

    ## Trick for dealing with singular MMt due to colinearity
    MMt <- MMt/max(MMt) + diag(0.95, nrow(MMt)) 
    #n <- nrow(MMt)
    #MMt<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*MMt)*MMt
    return(MMt)
  }

  .calcVC <- function(trait, currentX, MMt, ngpu)
  {
   cat(" performing emma.REMLE  ... \n")
    ## perform likelihood ratio test for variance component Var_g
    #res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100,ngpu=ngpu)
    res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, ngpu=ngpu)
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
   cat("\n\n\n                           FINAL MODEL  \n")
   cat(" ------------------------------------------------------------------------------------------  \n")
 }

.print_final  <- function(selected_loci, map, herit, extBIC )
{
  if (length(selected_loci) == 1 & any(is.na(selected_loci)))
  {
      cat("No significant marker-trait associations have been found. \n\n")
  }  else {
     .print_results(selected_loci=selected_loci, map=map, h=herit, extBIC=extBIC)
     cat(" ------------------------------------------------------------------------------------------  \n")
          cat("\n\n")
  }   ## end if else


}  ## end function print.finals

 .print_results <- function(itnum=NULL, selected_loci, map,h, extBIC)
 {  if(!is.null(itnum)){ 
       cat(" Significant marker-trait association found ... \n")
       cat(" Results after iteration ", itnum, "\n")
    }
    cat(sprintf("%15s  %10s        %10s     %10s        %10s \n", 
                 "SNP", "Chrm", "Map Pos",  "Col Number",       "extBIC"))
    for(ii in 1:length(selected_loci)){
       if(is.na(selected_loci[ii])){
       cat(sprintf("%15s  %10s        %10s        %8s           %-8.2f \n", 
        "Null Model", " ", " ", " ", extBIC[ii] ))
       }  else {
       cat(sprintf("%15s  %10s        %10s       %8s            %-8.2f \n", 
        map[[1]][selected_loci[ii]], map[[2]][selected_loci[ii]], as.character(map[[3]][selected_loci[ii]]), 
             selected_loci[ii], 
        extBIC[ii] ))
     }  ## end if else 
   }
    cat("\n")
 }



  .find_qtl <- function(geno, availmemGb, indxNA, selected_loci, MMt, invMMt, best_ve, best_vg, 
                       currentX, error_checking, ncpu, verbose, trait, ngpu )
  {
    ##  internal function: use only with AM
    if(verbose) cat(" Calculating H matrix   \n")
    H <- calculateH(MMt=MMt, varE=best_ve, varG=best_vg ) 
    if(verbose) cat(" Calculating P matrix - NOT GPU. \n")
    P <- calculateP(H=H, X=currentX ) 
    rm(H)
    gc()
 
    if (verbose)
              cat(" Calculating  square root of M %*% t(M) and it's inverse. \n")
    MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking, 
                              ngpu=ngpu ) 

    if (verbose)
            cat(" Calculating BLUPs for dimension reduced model. \n")
    hat_a <- calculate_reduced_a(varG=best_vg, P=P, 
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                       y=trait, verbose = verbose )   
     rm(P)
     gc()


    if (verbose) 
             cat(" Calculating variance of BLUPs for dimension reduced model. \n")
    var_hat_a    <- calculate_reduced_vara(X=currentX, varE=best_ve, varG=best_vg, invMMt=invMMt, 
                                                MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                                                verbose = verbose ) 

   
    if (verbose)
         cat(" Calculating BLUPs and their variances for full model. \n")
    ## not enough memory for this ......
     gc()
     ## load("everything.RData")   ## just for testing ... 
     a_and_vara  <- calculate_a_and_vara(maxmemGb=availmemGb, 
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
#' @description \code{AM} is used for multiple locus association mapping. It can simultaneously 
#' account for population stratification, familial relatedness, and nuisance fixed effects while 
#' detecting and mapping multiple marker-trait associations. It doesn't require any parameters to be tuned,
#' as with regularization technniques nor does it require a significance level or threshold to be set 
#' Those snp in strongest association with a trait are reported in a table of results. These 
#' snp map to the genomic regions containing the genes that are influencing the trait. 
#' @param trait  the name of the column in the phenotypic data file that contains the trait data. The name is case sensitive and must match exactly the column name in the phenotypic data file. 
#' @param feffects   a character vector containing the column names of 
#'                        the explanatory/independent variables in the phenotypic data file. If
#'                        not specified, only an overall mean will be fitted.
#' @param availmemGb a numeric value. It specifies the amount of available memory (in Gigabytes). 
#' This should be set to the maximum practical value of available memory for the analysis. 
#' @param geno   the R  object obtained from running \code{\link{ReadMarker}}. This must be specified. 
#' @param pheno  the R  object  obtained  from running \code{\link{ReadPheno}}. This must be specified.
#' @param map   the R object obtained from running \code{\link{ReadMap}}. If not specifed, a generic map will 
#'              be assumed. 
#' @param ncpu a numeric value for the number of CPU that are available for distributed computing.  The default is to determine the number of CPU automatically. 
#' @param ngpu   a integer value for the number of GPU available for computation.  The default
#'               is to assume there are no gpu available. 
#' @param  verbose      a logical value. When \code{TRUE}, extra output is returned 
#'  to the screen for monitoring progress. 
#' @param maxit     an integer value for the maximum number of forward steps to be performed. That is, it is the maximum number of 
#' genomic locations of interest that are to be identified. 
#' @param  error_checking a logical value. When \code{TRUE}, 
#' the numericial stability of the dimension reduction is checked. That is, individuals 
#' with near identical marker genotypes can cause numerical issues, and are reported 
#' when \code{error_checking} has been set to \code{TRUE}. 
#' @details
#'
#' \subsection{How to perform a basic AM analysis}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data is contained in the file "geno.txt" which is a plain space separated
#' text file with no column headings. The file is located in the current working directory or
#' the default directory for your R session.It contains numeric genotype values 0,1, and 2 for 
#' AA, AB, and BB, respectively}
#' \item{}{the phenotypic data is contained in the file "pheno.txt" which is a plain space
#' separated text file with a single trait and no explanatory variables. The file has the 
#' column heading "y". The file is located in the current working directory.}
#' \item{}{there is no map data}
#' }
#'
#'  To analyse these data, we would run the following:
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename="geno.txt", AA=0, AB=1, BB=2, type="text")
#'   
#'   pheno_obj <- ReadPheno(filename="pheno.txt")
#'
#'   res <- AM(trait="y", geno=geno_obj, pheno=pheno_obj)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{How to perform a more complicated AM analysis}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data is contained in the file "geno.ped" which is a PLINK ped file. See
#' \code{\link{ReadMarker}} for details. The file is located in /my/dir. Lets assume 
#' the file is large but our machine has 32Gbytes of RAM.}
#' \item{}{the phenotypic data is contained in the file "pheno.txt" which is a plain space
#' separated text file with  six columns. The first column is a trait and is labelled "y1".
#' The second column is another trait and is laballed "y2". The third and fourth columns 
#' are nuisance variables and are labelled "cov1" and "cov2". The fifth and sixth columns
#' are the first two principal components to account for population substructure and are 
#' labelled "pc1" and "pc2".
#' The file is located in /my/dir}
#' \item{}{the map data is contained in the file "map.txt" and is also located in 
#'  /my/dir}
#' \item{}{An AM analysis is performed where the trait of interest is "trait2", 
#' the fixed effects to be included in the analysis are "cov1", "cov2", "pc1", and "pc2", 
#' and the available memory is to be set to 32Gbytes.}
#' } 
#'
#'  To analyse these data, we would run the following:
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename="/my/dir/geno.ped", type="PLINK", availmemGb=32)
#'   
#'   pheno_obj <- ReadPheno(filename="/my/dir/pheno.txt")
#'
#'   map_obj   <- ReadMap(filename="/my/dir/map.txt")
#'
#'   res <- AM(trait="trait2", feffects=c("cov1", "cov2", "pc1", "pc2"), 
#'             geno=geno_obj, pheno=pheno_obj, map=map_obj, availmemGb=32)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{Dealing with Missing Values}{
#'  
#'
#' }



#' The genotypic file must not contain missing genotypes.
#'
#' The phenotypic file can contain missing information, coded as \code{NA}, but only 
#' in the trait data. The fixed effects data (i.e. the explanatory variables)
#' cannot contain any missing data.  The phenotypic file is allowed to contain multiple traits
#' and fixed effects.  
#'  
#' The trait data and fixed effects data are  specified by setting \code{trait} 
#'  and \code{feffects}, respectively. \code{trait} can only contain a single
#'  character string for the trait name.  \code{feffects} can be a character vector 
#'  if multiple fixed effects are to be included in the linear mixed model. Whether the 
#'  fixed effects are treated as a covariate or factor is dependent upon the class of 
#'  the associated data columns in the data frame obtained from \code{\link{ReadPheno}}. 
#'
#'    STILL BEING WRITTEN .... 
#' The \code{AM} function is an R/Rcpp implementation of multi-locus whole-genome 
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
#' @seealso \code{\link{ReadMarker}}, \code{\link{ReadPheno}}, and \code{\link{ReadMap}}
#'
#' @return
#' something here .... 
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
#'                            feffects = c("cov1", "cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
#'
#'
#'
#'
#'
AM <- function(trait=NULL, 
               feffects  = NULL,
               availmemGb=8, 
               geno=NULL, 
               pheno=NULL, 
               map = NULL,
               ncpu=detectCores(),
               ngpu=0,
               verbose=FALSE,
               maxit=20,
               error_checking=FALSE 
               ){

 ## Core function for performing whole genome association mapping with EMMA
 ## Args
 ## ncpu        number of cores available for computation
 ## memoryGb        maximum amount of working memory available for computation
 ## pheno           data frame 
 ##                 remaining columns are explanatory variables to include in the model. If a numeric vector, then it 
 ##                 is only a response to be fitted. 
 ## geno            if geno is a matrix or data frame, then the user has not ReadMarker and a bin packed file
 ##                 has not been created. If it is a character string, then it is the file location of the binary packed files. 
 ## error_checking  when true, it performs some checks of the calculations
 ## maxit           maximum number of qtl to include in the model
 ## ngpu            number of GPU available for computation

 ## check parameter inputs
 check.inputs.mlam(ncpu, availmemGb, trait, feffects, 
                     map, pheno, geno )
 ## ADD CHECK TO MAKE SURE PHENO AND GENO NROWS ARE CORRECT - NOT BEING CHECKED AT THE MOMENT

 selected_loci <- NA
 new_selected_locus <- NA
 herit <- vector("numeric", 0)
 extBIC <- vector("numeric", 0)
 ## assign trait 
 trait <-  pheno[[trait]]
 
 ## check for NA's in explanatory variables. 
 ## If any, set individual's trait value to NA
 ## This means this individual will later be removed. 
 if(!is.null(feffects)){
  mat.of.NA  <- which(is.na(pheno[, feffects]), arr.ind=TRUE)
  if(dim(mat.of.NA)[1]>0){
    trait[unique(mat.of.NA[,1])] <- NA
  }
 }

 ## check for NA's in trait
 indxNA <- check.for.NA.in.trait(trait=trait)



  if(ngpu > 0 ){
# library(rcppMagmaSYEVD)
## caters for the two ways data can be inputed into AMplus
     cat(" oooooooooooooooo \n")
       rcppMagmaSYEVD::RunServer( matrixMaxDimension=geno[["dim_of_bin_M"]][1],  numGPUsWanted=ngpu, memName="/syevd_mem", semName="/syevd_sem", print=0)
  }


 ## remove missing observations from trait
 if(length(indxNA)>0){
    trait <- trait[-indxNA]
 }

 ## build design matrix currentX
 cat(" Forming currentX \n")
 Args <- list(pheno=pheno, geno=geno, indxNA=indxNA, feffects=feffects, verbose=verbose )
currentX <- do.call(.build_design_matrix, Args) 

 ## print tile
 .print_title()

 ## Initialization
 continue <- TRUE
 itnum <- 1


 while(continue){
   cat(" Performing iteration ... ", itnum, "\n")

   ## based on selected_locus, form model matrix X
  currentX <- constructX(currentX=currentX, loci_indx=new_selected_locus,
                          dim_of_bin_M=geno[["dim_of_bin_M"]],
                          indxNA = indxNA,
                          map=map, availmemGb = availmemGb)  



    ## calculate Ve and Vg
    Args <- list(geno=geno,availmemGb=availmemGb,
                    ncpu=ncpu,selected_loci=selected_loci,
                    verbose=verbose, indxNA=indxNA)

    if(itnum==1){
         cat(" Calculate MMt \n")   
         MMt <- do.call(.calcMMt, Args)  

        invMMt <- chol2inv(chol(MMt))   ## doesn't use GPU
        gc()
    } 
      vc <- .calcVC(trait=trait, currentX=currentX,MMt=MMt, ngpu=ngpu) 
    gc()
    best_ve <- vc[["ve"]]
    best_vg <- vc[["vg"]]


    ## Calculate extBIC
    cat(" Calculate extBIC ... \n")
    new_extBIC <- .calc_extBIC(trait, currentX,MMt, geno, verbose) 
    h <- best_vg/(best_vg + best_ve)
    gc()

    ## set vectors herit and extBIC
    herit <- c(herit, h)
    extBIC <- c(extBIC, new_extBIC)


    ## Print findings to screen

   .print_results(itnum, selected_loci, map, herit, extBIC)

   ## increased June 1, 2016 because it was stopping short. 
   ## Select new locus if h > 0.01 | or extBIC is still decreasing 
   ## if(h > 0.01 | which(extBIC==min(extBIC))==length(extBIC) ){
   if(which(extBIC==min(extBIC))==length(extBIC) ){  ## new way of stoppint based on extBIC only
     ## find QTL
     ARgs <- list(geno=geno,availmemGb=availmemGb, indxNA=indxNA, selected_loci=selected_loci,
                 MMt=MMt, invMMt=invMMt, best_ve=best_ve, best_vg=best_vg, currentX=currentX,
                 error_checking=error_checking,
                 ncpu=ncpu, verbose=verbose, trait=trait, ngpu=ngpu)
      cat(" finding QTL ..,. \n")
      new_selected_locus <- do.call(.find_qtl, ARgs)  ## memory blowing up here !!!! 
     gc()
     selected_loci <- c(selected_loci, new_selected_locus)

   }  else {
     cat( " while terminating loop ... \n")
     ## terminate while loop, h near 0
     continue <- FALSE
     .print_header()
     .print_final(selected_loci, map, herit, extBIC)
     sigres <- .form_results(trait, selected_loci, map,  feffects, 
                     indxNA, ncpu, availmemGb, verbose, herit, extBIC )   
   }  ## end if else


   itnum <- itnum + 1
   ## alternate stopping rule - if maxit has been exceeded.
    if(itnum > maxit){
         continue <- FALSE 
         .print_header()
         ## need to remove the last selected locus since we don't go on and calucate its H and extBIC 
         ## under this new model. 
         .print_final(selected_loci[-length(selected_loci)], map, herit, extBIC)
         sigres <- .form_results(trait, selected_loci[-length(selected_loci)], map,  feffects, 
                     indxNA, ncpu, availmemGb, verbose, herit, extBIC )   
    }
 
  }  ## end while continue





return( sigres )

} ## end AM















