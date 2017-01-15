# This software is distributed under the GNU General Public License.
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 



doquiet <- function(dat, num_markers, lab){
     ## a diagnostic function for printing the contents of matrix or vector
     ## used for error checking

     if(dim(dat)[1] == 1 || dim(dat)[2] ==1 )
         dat <- as.numeric(dat)

     if(class(dat)=="matrix"){
          ### if dat is a matrix

         if(num_markers > 0){
           cat(" Dimension of ", lab, " is", dim(dat), "\n")
           cat(" First few rows and ", num_markers, "columns of ", lab, " are: \n")
           if(nrow(dat) > 5 && ncol(dat) > num_markers)
               print(dat[1:5, 1:num_markers])
           if(nrow(dat) <=5  &&  ncol(dat) > num_markers)
               print(dat[1:nrow(dat), 1:num_markers])
           if(nrow(dat) > 5  &&  ncol(dat) <=  num_markers)
               print(dat[1:5, 1:ncol(dat)])
           if(nrow(dat) <= 5  &&  ncol(dat) <=  num_markers)
               print(dat[1:nrow(dat), 1:ncol(dat)])
           cat("\n\n")
         }
     } ## end if class(dat)

     if(class(dat)=="numeric" || class(dat)=="vector"){
       if(num_markers > 0){
          cat(" Length of ", lab, "is", length(dat), "\n")
          cat(" The first ", num_markers, "elements of the vector are ", lab, "\n")
          if(length(dat) > num_markers)
             print(dat[1:num_markers])
          if(length(dat) <= num_markers)
             print(dat[1:length(dat)])
       cat("\n\n")
       }
    }

    if(!(class(dat)=="matrix" || class(dat)=="vector" || class(dat)=="numeric"))
      cat(" Internal error in doquiet. dat not matrix or vector or numeric. \n")

}

.form_results <- function(trait, selected_loci, map,  feffects, indxNA,
                           ncpu, availmemGb, quiet,  extBIC )
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
                    quiet=quiet,
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
                    quiet=quiet,
                    extBIC=extBIC)
  }
return(sigres)
}

.print_title <- function(){
    ## internal fuction: use only in AM function
    ## title
    cat("\n\n\n\n")
cat("                    Multiple-Locus Association Mapping\n")
cat("                            Version 1.0 \n\n")
cat(" \n")
cat("   . ,-\"-.   ,-\"-. ,-\"-.   ,-\"-. ,-\"-. ,-\"-. ,-\"-.   ,-\"-. ,-\"-.    \n")  
cat("    X | | \\ / | | X | | \\ / | | X | | \\ / | | X | | \\ / | | X | | \\ /   \n")
cat("   / \\| | |X| | |/ \\| | |X| | |/ \\| | |X| | |/ \\| | |X| | |/ \\| | |X|   \n")
cat("      `-!-' `-!-\"   `-!-' `-!-'   `-!-' `-!-\"   `-!-' `-!-'   `-!-' `-     \n\n\n")

}


.build_design_matrix <- function(pheno=NULL, geno=NULL, indxNA=NULL, feffects=NULL, quiet=0  )
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

 if (quiet > 0){
   cat("Dimension of design matrix, before addition of marker fixed effects is ", nrow(Xmat), "rows and ", ncol(Xmat), "columns.\n") 
 }

if(!is.matrix(Xmat))
   Xmat <- matrix(data=Xmat, ncol=1)

  return(Xmat)
}


.calcMMt <- function(geno, availmemGb, ncpu, selected_loci, quiet, indxNA)
  {
    ## internal function: used only in multilocus_loci_am and summaryam
    ## values passed by environments
    MMt <- calculateMMt(geno=geno[["asciifileM"]], availmemGb=availmemGb, 
                           ncpu=ncpu, 
                           dim_of_ascii_M = geno[["dim_of_ascii_M"]], 
                           selected_loci=selected_loci, quiet = quiet) 
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
    ## perform likelihood ratio test for variance component Var_g
    #res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100,ngpu=ngpu)
    res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, ngpu=ngpu)
    return(list("vg"=res_full$vg, "ve"=res_full$ve))

  }

 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, quiet=FALSE)
 { 
   ## smallest extBIC and BIC is best
   ## internal function: use in multiple_loucs_am only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100)
   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

   extBIC <- BIC + 2 * lchoose(geno$dim_of_ascii_M[2], ncol(currentX) - 1)  

    return(extBIC)
 }



 .print_header <- function(){
   cat("\n\n\n                           Final  Results  \n")
   cat(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
 }

.print_final  <- function(selected_loci, map,  extBIC )
{
  if (length(selected_loci) == 1 & any(is.na(selected_loci)))
  {
      cat("No significant marker-trait associations have been found. \n\n")
  }  else {
     .print_results(selected_loci=selected_loci, map=map,  extBIC=extBIC)
          cat("\n\n")
  }   ## end if else


}  ## end function print.finals

 .print_results <- function(itnum=NULL, selected_loci, map, extBIC)
 {  if(!is.null(itnum)){ 
       cat(" Significant marker-trait association found. \n\n")
       cat(" New results after iteration ", itnum, "are \n\n")
    }
    cat(sprintf("%15s  %10s        %10s     %10s        %10s \n", 
                 "SNP", "Chrm", "Map Pos",  "Col Number",       "extBIC"))
    cat(sprintf("%15s  %10s        %10s     %10s        %10s \n", 
                 "-----", "------", "---------",  "-----------",       "---------"))

    for(ii in 1:length(selected_loci)){
       if(is.na(selected_loci[ii])){
       cat(sprintf("%15s  %10s        %10s        %8s           %-8.2f \n", 
        "Null Model", " ", " ", " ", extBIC[ii] ))
       }  else {
       cat(sprintf("%15s  %10s        %10s       %8s            %-8.2f \n", 
        map[[1]][selected_loci[ii]], map[[2]][selected_loci[ii]], as.character(map[[3]][selected_loci[ii]]), 
             selected_loci[ii], extBIC[ii] ))
     }  ## end if else 
   }
    cat("\n\n\n\n")
 }



  .find_qtl <- function(geno, availmemGb, indxNA, selected_loci, MMt, invMMt, best_ve, best_vg, 
                       currentX,  ncpu, quiet, trait, ngpu )
  {
    ##  internal function: use only with AM
    if(quiet > 0){
       cat(" quiet =", quiet, ": beginning calculation of H matrix. \n")
    }
    H <- calculateH(MMt=MMt, varE=best_ve, varG=best_vg ) 
    doquiet(dat=H, num_markers=quiet, lab="H")

    if(quiet>0){
       cat(" quiet =", quiet, ": beginning calculation of P matrix. \n")
    }
    P <- calculateP(H=H, X=currentX ) 
    doquiet(dat=P, num_markers=quiet, lab="P")
    rm(H)
    gc()
 
    
    if(quiet > 0){
      cat(" quiet = ", quiet, ": beginning calculation of the square root of MMt and its inverse. \n")
    }
    ## artifact from old code but kept it anyway. Looks at the stability of the MMt calculation 
    ## especially if there are near identical rows of data in M
    error_checking <- FALSE
    if (quiet > 0)
       error_checking <- TRUE
    MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking, 
                              ngpu=ngpu ) 

    doquiet(dat=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], num_markers=quiet, lab="sqrt(M %*% M^t)")
    doquiet(dat=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]], num_markers=quiet, lab="sqrt(M %*% M^t)^-1")

    if(quiet > 0){
      cat(" quiet =", quiet, ": beginning calculation of the BLUP estimates for dimension reduced model. \n")
    }
    hat_a <- calculate_reduced_a(varG=best_vg, P=P, 
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                       y=trait, quiet = quiet )   
    doquiet(dat=hat_a, num_markers=quiet, lab="BLUPs")


     rm(P)
     gc()

    if(quiet > 0){
      cat(" quiet = ", quiet, ": beginning calculation of the standard errors  of BLUP estimates for dimension reduced model. \n")
    }

    var_hat_a    <- calculate_reduced_vara(X=currentX, varE=best_ve, varG=best_vg, invMMt=invMMt, 
                                                MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                                                quiet = quiet ) 
    doquiet(dat=var_hat_a, num_markers=quiet, lab="SE of BLUPs")


   
     gc()
     ## load("everything.RData")   ## just for testing ... 
    if(quiet > 0){
      cat(" quiet = ", quiet, ": beginning calculation of BLUPS and their standard errors for full model. \n")
    }

     a_and_vara  <- calculate_a_and_vara(maxmemGb=availmemGb, 
                                            dims=geno[["dim_of_ascii_M"]],
                                            selectedloci = selected_loci,
                                            invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]],
                                            transformed_a=hat_a, 
                                            transformed_vara=var_hat_a,
                                            indxNA = indxNA,
                                            quiet=quiet) 

     doquiet(dat=a_and_vara[["a"]], num_markers=quiet, lab="BLUPs for full model")
     doquiet(dat=a_and_vara[["vara"]], num_markers=quiet, lab="SE of BLUPs for full model")


  
    ## outlier test statistic
    if (quiet > 0) 
        cat(" quiet = ", quiet, ": beginning calculation of  outlier test statistics. \n")
    tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
    doquiet(dat=tsq, num_markers=quiet, lab="outlier test statistic")

    indx <- which(tsq == max(tsq, na.rm=TRUE))   ## index of largest test statistic. However, need to account for other loci 
                                         ## already having been removed from M which affects the indexing

    ## taking first found qtl
    indx <- indx[1]

    orig_indx <- seq(1, geno[["dim_of_ascii_M"]][2])  ## 1:ncols
    return(orig_indx[indx])
}

#' @title multiple-locus association mapping 
#' @description \code{AM} is used for multiple locus association mapping. It can simultaneously 
#' account for population stratification, familial relatedness, and nuisance fixed effects while 
#' detecting and mapping multiple marker-trait associations. It doesn't require any parameters to be tuned,
#' as with regularization technniques nor does it require a significance level or threshold to be set. 
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
#' @param  quiet      an integer value specifying the number of marker loci for which diagnostic information is 
#' to be printed to the screen. This is useful for error checking. 
#' @param maxit     an integer value for the maximum number of forward steps to be performed. That is, it is the maximum number of 
#' genomic locations of interest that are to be identified. 
#' @details
#'
#' \subsection{How to perform a basic AM analysis}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data is contained in the file "geno.txt" which is a plain space separated
#' text file with no column headings. The file is located in the current working directory. 
#' It contains numeric genotype values 0, 1, and 2 for 
#' AA, AB, and BB, respectively}
#' \item{}{the phenotypic data is contained in the file "pheno.txt" which is a plain space
#' separated text file with a single trait and no explanatory variables. The file has the 
#' column heading "y". The file is located in the current working directory.}
#' \item{}{there is no map data}
#' }
#'
#'  To analyse these data, we would use the following functions:
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
#' the file is large but our machine has 32 Gbytes of RAM.}
#' \item{}{the phenotypic data is contained in the file "pheno.txt" which is a plain space
#' separated text file with  six columns. The first column is a trait and is labelled "y1".
#' The second column is another trait and is laballed "y2". The third and fourth columns 
#' are nuisance variables and are labelled "cov1" and "cov2". The fifth and sixth columns
#' are the first two principal components to account for population substructure and are 
#' labelled "pc1" and "pc2".
#' The file is located in /my/dir.}
#' \item{}{the map data is contained in the file "map.txt", is also located in 
#'  /my/dir, and the first row has the column headings.}
#' \item{}{An AM analysis is performed where the trait of interest is "y2", 
#' the fixed effects to be included in the analysis are "cov1", "cov2", "pc1", and "pc2", 
#' and the available memory is to be set to 32 Gbytes.}
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
#'   res <- AM(trait="y2", feffects=c("cov1", "cov2", "pc1", "pc2"), 
#'             geno=geno_obj, pheno=pheno_obj, map=map_obj, availmemGb=32)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{Dealing with missing marker data}{
#'
#' In running \code{\link{AM}}, we assume there are no missing marker genotypes. 
#' In fact, \code{ReadMarker} will generate errors if there are missing marker data.  
#' However, most genome-wide data sets have missing genotypes (sometimes by design).  
#' Ideally, a genotype imutation program such as  BEAGLE, MACH, fastPHASE, or PHASE2, should be 
#' used to impute the missing marker data. 
#'
#'
#' Alternately, for quick results, remove those loci with a high proportion of missing data (> 10\%) 
#' and  replace the remaining missing marker genotypes with  heterozygote genotypes. 
#' Since we assume an additive model in \code{AM}, this will not cause false positives but it 
#' can reduce power. We found though that the loss in  power is minimal if the 
#' proportion of missing data is low.  See George and Cavanagh (2015) for details.  
#' }
#'
#' \subsection{Dealing with missing trait data}{
#'
#'  \code{AM} deals automatically with individuals with missing trait data. 
#' These individuals are removed  from the analysis and a warning message is generated.
#' }
#' 
#' \subsection{Dealing with missing explanatory variable values}{
#'
#' \code{AM} deals automatically with individuals with missing explanatory variable values. 
#' These individuals are removed from the analysis and a warning message is generated
#' }
#'
#' \subsection{Error Checking}{
#'
#' \code{quiet} specifies the number of marker loci for which diagnostic information 
#' is to be printed to the screen.  When \code{quiet} is non-zero, the contents of important matrices and 
#' vectors are printed. This can be useful for diagnosing problems with the 
#' input data, especially when you compare the contents of the matrices/vectors for data that loads correctly 
#' to data that is causing errors. 
#'
#' Setting \code{quiet} to an integer value can also be useful for monitoring progress when analysing large 
#' data sets.   
#'}
#'
#'
#'
#' @references George AW and Cavanagh C. 2015. Genome-wide Association Mapping in Plants. 
#' Theorectical and Applied Genetics 128: 1163-1174.
#'
#'
#'
#' @seealso \code{\link{ReadMarker}}, \code{\link{ReadPheno}}, and \code{\link{ReadMap}}
#'
#' @return
#' A list with the following components:
#' \describe{
#'\item{trait}{column name of the trait being used by AM}
#'\item{feffects}{column names of the explanatory variables being used by AM}
#'\item{indxNA}{numeric vector containing the line numbers of those individuals with missing phenotypic data for the 
#' trait and explanatory variables being used by AM}
#' \item{Mrk}{character vector with the marker names of those loci found to be in significant association with the trait}
#' \item{Chr}{character vector with the chromosome names for the loci found to be in significant association with the trait}
#' \item{Pos}{numeric vector with the map positions of the loci found to be in  significant association with the trait}
#' \item{Indx}{column number in the marker file of the loci found to be in  significant association with the trait}
#' \item{ncpu}{number of cpu used for the calculations}
#' \item{availmemGb}{amount of RAM in Gbytes that has been set by the user}
#' \item{quiet}{the number of markers for which diagnostic information is to be printed.}
#' \item{extBIC}{numeric vector with the extended BIC values for the loci  found to be in  significant association with the trait}
#'}
#'
#' @examples
#'   #-------------------------
#'   #  Example  
#'   #------------------------
#'
#'   # read the map 
#'   #~~~~~~~~~~~~~~
#'   
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file("extdata", "map.txt", 
#'                                    package="AMplus")
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'  # to look at the first few rows of the map file
#'  head(map_obj)
#'
#'   # read marker data
#'   #~~~~~~~~~~~~~~~~~~~~
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8Gbytes
#'   complete.name <- system.file("extdata", "geno.ped", 
#'                                      package="AMplus")
#'   geno_obj <- ReadMarker(filename=complete.name,  type="PLINK", availmemGb=8) 
#'  
#'   # read phenotypic data
#'   #~~~~~~~~~~~~~~~~~~~~~~~
#'
#'   # Read in a plain text file with data on a single trait and two covariates
#'   # The first row of the text file contains the column names "trait", "cov1", and "cov2". 
#'   complete.name <- system.file("extdata", "pheno.txt", package="AMplus")
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'   # Perform multiple-locus genome-wide association mapping 
#'   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'   res <- AM(trait = "trait",
#'                            feffects = c("cov1", "cov2"),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
#'
#'  # Performing multiple-locus genome-wide association mapping with a model 
#'  #    with no fixed effects except for an intercept. 
#'  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'  
#'   res <- AM(trait = "trait",
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
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
               quiet=0,
               maxit=20
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
 ## maxit           maximum number of qtl to include in the model
 ## ngpu            number of GPU available for computation

 ## check parameter inputs

 ## print tile
 .print_title()


 error.code <- check.inputs.mlam(ncpu=ncpu , availmemGb=availmemGb, colname.trait=trait, colname.feffects=feffects, 
                     map=map, pheno=pheno, geno=geno )
 if(error.code)
    stop("AM has terminted with errors.", call. = FALSE)




 ## checking if map is present. If not, generate a fake map. 
 if(is.null(map)){
   if(quiet > 0){
     cat(" Map file has not been supplied. An artifical map is being created but this map is not used in the analysis. \n")
     cat(" It is only used for the reporting of results. \n")
   }
   ## map has not been supplied. Create own map
   map <- data.frame(SNP=paste("M", 1:geno[["dim_of_ascii_M"]][2], sep=""), 
                     Chr=rep(1, geno[["dim_of_ascii_M"]][2]), 
                     Pos=1:geno[["dim_of_ascii_M"]][2])
  }

 ## check that the number of rows in the map file match the number of columns in the geno file
 if (geno[["dim_of_ascii_M"]][2] != nrow(map)){
   cat(" Error: There is a differing number of loci read in by ReadMarker and ReadMap functions. \n")
   cat("         The number of marker loci read in by ReadMarker() is ", geno[["dim_of_ascii_M"]][2], "\n")
   cat("        The number of marker loci in  the marker map is  ", nrow(map), "\n") 
   stop(" AM has terminatated with errors.", call. = FALSE)
 }


 ## check that the number of rows in the phenotype file match the number of rows in the geno file
 if (geno[["dim_of_ascii_M"]][1] != nrow(pheno)){
   cat(" Error: There is a differing number  of rows read in by ReadMarker and ReadPheno functions. \n")
   cat("         The number of rows read in by ReadMarker() is ", geno[["dim_of_ascii_M"]][1], "\n")
   cat("        The number of rows  read in by ReadPheno is  ", nrow(map), "\n") 
   stop(" AM has terminatated with errors.", call. = FALSE)
 }




 selected_loci <- NA
 new_selected_locus <- NA
 extBIC <- vector("numeric", 0)
 ## assign trait 
 trait <-  pheno[[trait]]
 
 ## check for NA's in explanatory variables. 
 ## If any, set individual's trait value to NA
 ## This means this individual will later be removed. 
 if(!is.null(feffects)){
  mat.of.NA  <- which(is.na(pheno[, feffects]), arr.ind=TRUE)
  if(!is.null(dim(mat.of.NA)[1]) ){
     if(dim(mat.of.NA)[1]>0){
       trait[unique(mat.of.NA[,1])] <- NA
     }
  }
 }

 ## check for NA's in trait
 indxNA <- check.for.NA.in.trait(trait=trait)


  ## setting up gpu sevrver
  if(ngpu > 0 ){
     if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
        library(rcppMagmaSYEVD)
         rcppMagmaSYEVD::RunServer( matrixMaxDimension=geno[["dim_of_ascii_M"]][1],  numGPUsWanted=ngpu, memName="/syevd_mem", semName="/syevd_sem", print=0)
     } 
  }


 ## remove missing observations from trait
 if(length(indxNA)>0){
    trait <- trait[-indxNA]

    if(quiet > 0){
     cat(" The following rows are being removed from pheno due to missing data: \n")
     cat("             ", indxNA, "\n\n")
    }

 }


 ## build design matrix currentX
currentX <- .build_design_matrix(pheno=pheno, geno=geno, indxNA=indxNA, feffects=feffects, quiet=quiet )

 ## Initialization
 continue <- TRUE
 itnum <- 1


 while(continue){
  cat("\n\n Iteration" , itnum, ": Searching for most significant marker-trait association\n\n")
   ## based on selected_locus, form model matrix X
  currentX <- constructX(currentX=currentX, loci_indx=new_selected_locus,
                          dim_of_ascii_M=geno[["dim_of_ascii_M"]],
                          indxNA = indxNA,
                          map=map, availmemGb = availmemGb)  



    ## calculate Ve and Vg
    Args <- list(geno=geno,availmemGb=availmemGb,
                    ncpu=ncpu,selected_loci=selected_loci,
                    quiet=quiet, indxNA=indxNA)

    if(itnum==1){
        if(quiet>0)
           cat(" quiet=FALSE: calculating M %*% M^t. \n")
         MMt <- do.call(.calcMMt, Args)  



         doquiet(dat=MMt, num_markers=quiet, lab="M%*%M^t")
        invMMt <- chol2inv(chol(MMt))   ## doesn't use GPU
        gc()
    } 
    if(quiet>0){
      cat(" Calculating variance components for multiple-locus model. \n")
    }
    vc <- .calcVC(trait=trait, currentX=currentX,MMt=MMt, ngpu=ngpu) 
    gc()
    best_ve <- vc[["ve"]]
    best_vg <- vc[["vg"]]

   if(quiet>0){
      cat(" Residual variance estimate is ", best_ve, "\n")
      cat(" Polygenic variance estimate is ", best_vg, "\n")
   }


    ## Calculate extBIC
    new_extBIC <- .calc_extBIC(trait, currentX,MMt, geno, quiet) 
    gc()

    ## set vector extBIC
    extBIC <- c(extBIC, new_extBIC)


    ## Print findings to screen
   .print_results(itnum, selected_loci, map,  extBIC)
   ## Select new locus if extBIC is still decreasing 
   if(which(extBIC==min(extBIC))==length(extBIC) ){  ## new way of stoppint based on extBIC only
     ## find QTL
     ARgs <- list(geno=geno,availmemGb=availmemGb, indxNA=indxNA, selected_loci=selected_loci,
                 MMt=MMt, invMMt=invMMt, best_ve=best_ve, best_vg=best_vg, currentX=currentX,
                 ncpu=ncpu, quiet=quiet, trait=trait, ngpu=ngpu)
      new_selected_locus <- do.call(.find_qtl, ARgs)  ## memory blowing up here !!!! 
     gc()
     selected_loci <- c(selected_loci, new_selected_locus)

   }  else {
     ## terminate while loop, 
     continue <- FALSE
#     .print_header()
#     .print_final(selected_loci, map, extBIC)
#     sigres <- .form_results(trait, selected_loci, map,  feffects, 
#                     indxNA, ncpu, availmemGb, quiet,  extBIC )   
   }  ## end if else


   itnum <- itnum + 1
   ## alternate stopping rule - if maxit has been exceeded.
    if(itnum > maxit){
         continue <- FALSE 
         .print_header()
         ## need to remove the last selected locus since we don't go on and calucate its H and extBIC 
         ## under this new model. 
         .print_final(selected_loci[-length(selected_loci)], map, extBIC)
         sigres <- .form_results(trait, selected_loci[-length(selected_loci)], map,  feffects, 
                     indxNA, ncpu, availmemGb, quiet,  extBIC )   
    }
 
  }  ## end while continue

if( itnum > maxit){
    .print_header()
    .print_final(selected_loci, map,  extBIC)
    sigres <- .form_results(trait, selected_loci, map,  feffects, 
                     indxNA, ncpu, availmemGb, quiet,  extBIC )   

} else {
    ## remove last selected_loci as for this locus, the extBIC went up
    if(length(selected_loci)>1){
        .print_header()
        .print_final(selected_loci[-length(selected_loci)], 
                     map, 
                     extBIC[-length(selected_loci)])
        sigres <- .form_results(trait, selected_loci[-length(selected_loci)], map,  feffects, 
                         indxNA, ncpu, availmemGb, quiet, 
                         extBIC[-length(selected_loci)] )   
    } else {
        .print_header()
        .print_final(selected_loci, map, extBIC)
        sigres <- .form_results(trait, selected_loci, map,  feffects, 
                         indxNA, ncpu, availmemGb, quiet, extBIC )   
   }  ## end inner  if(lenth(selected_locus)>1)
}  ## end if( itnum > maxit)


return( sigres )

} ## end AM















