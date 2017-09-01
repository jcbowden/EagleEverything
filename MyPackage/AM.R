
#' @title multiple-locus association mapping 
#' @description \code{AM} performs  association mapping within a multiple-locus linear mixed model framework. 
#' \code{AM}  finds the best set of 
#' marker loci in strongest association with a trait while simultaneously accounting for any fixed effects and the genetic background.     
#' @param trait  the name of the column in the phenotype data file that contains the trait data. The name is case sensitive and must match exactly the column name in the phenotype data file. 
#' @param fformula   the right hand side formula for the fixed effects.   See below for details. 
#'                        If
#'                        not specified, only an overall mean will be fitted.
#' @param availmemGb a numeric value. It specifies the amount of available memory (in Gigabytes). 
#' This should be set to the maximum practical value of available memory for the analysis. 
#' @param geno   the R  object obtained from running \code{\link{ReadMarker}}. This must be specified. 
#' @param pheno  the R  object  obtained  from running \code{\link{ReadPheno}}. This must be specified.
#' @param map   the R object obtained from running \code{\link{ReadMap}}. If not specified, a generic map will 
#'              be assumed. 
#' @param ncpu a integer  value for the number of CPU that are available for distributed computing.  The default is to determine the number of CPU automatically. 
#' @param ngpu   a integer value for the number of gpu available for computation.  The default
#'               is to assume there are no gpu available.  This option has not yet been implemented.
#' @param  quiet      a logical value. If set to \code{TRUE}, additional runtime output is printed. 
#' This is useful for error checking and monitoring the progress of a large analysis. 
#' @param maxit     an integer value for the maximum number of forward steps to be performed.  This will rarely need adjusting. 
#' @details
#'
#' \subsection{How to perform a basic AM analysis}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data are contained in the file geno.txt which is a plain space separated
#' text file with no column headings. The file is located in the current working directory. 
#' It contains numeric genotype values 0, 1, and 2 for snp genotypes
#' AA, AB, and BB, respectively. It also contains the numeric value X for a missing genotype. }
#' \item{}{the phenotype data is contained in the file pheno.txt which is a plain space
#' separated text file containing a single column with the trait data. The first row of the file 
#' has the column heading 'y'. 
#' The file is located in the current working directory.}
#' \item{}{there is no map data.}
#' }
#'
#'  To analyse these data, we would use the following three functions:
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename='geno.txt', AA=0, AB=1, BB=2, type="text", missing='X')
#'   
#'   pheno_obj <- ReadPheno(filename='pheno.txt')
#'
#'   res <- AM(trait='y', geno=geno_obj, pheno=pheno_obj)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{How to perform a more complicated AM analysis}{
#'
#' Suppose, 
#' \itemize{
#' \item{}{the snp data are contained in the file geno.ped which is a 'PLINK' ped file. See
#' \code{\link{ReadMarker}} for details. The file is located in /my/dir. Let's assume 
#' the file is large, say 50 gigabytes,   and our computer only has 32 gigabytes of RAM.}
#' \item{}{the phenotype data is contained in the file pheno.txt which is a plain space
#' separated text file with  six columns. The first row of the file contains the column headings. 
#' The first column is a trait and is labeled y1.
#' The second column is another trait and is labeled y2. The third and fourth columns 
#' are nuisance variables and are labeled cov1 and cov2. The fifth and sixth columns
#' are the first two principal components to account for population substructure and are 
#' labeled pc1 and pc2. The file contains missing data that are coded as 99. 
#' The file is located in /my/dir.}
#' \item{}{the map data is contained in the file map.txt, is also located in 
#'  /my/dir, and the first row has the column headings.}
#' \item{}{An 'AM' analysis is performed where the trait of interest is y2, 
#' the fixed effects part of the model is cov1 + cov2 + pc1 + pc2, 
#' and the available memory is set to 32 gigabytes.}
#' } 
#'
#'  To analyse these data, we would run the following:
#' \preformatted{
#'   geno_obj <-  ReadMarker(filename='/my/dir/geno.ped', type='PLINK', availmemGb=32)
#'   
#'   pheno_obj <- ReadPheno(filename='/my/dir/pheno.txt', missing=99)
#'
#'   map_obj   <- ReadMap(filename='/my/dir/map.txt')
#'
#'   res <- AM(trait='y2', fformula=c('cov1 + cov2 + pc1 + pc2'), 
#'             geno=geno_obj, pheno=pheno_obj, map=map_obj, availmemGb=32)
#' }
#' A table of results is printed to the screen and saved in the R object \code{res}. 
#'}
#'
#' \subsection{Dealing with missing marker data}{
#'
#' \code{AM} can tolerate some missing marker data. However, ideally, 
#' a specialized genotype imputation program such as  'BEAGLE', 'MACH', 'fastPHASE', or 'PHASE2', should be 
#' used to impute the missing marker data before being read into 'Eagle'.  
#'
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
#' Most errors occur when reading in the data. However, as an extra precaution, if \code{quiet=TRUE}, then additional 
#' output is printed during the running of \code{AM}. If \code{AM} is failing, then this output can be useful for diagnosing 
#' the problem. 
#'}
#'
#'
#'
#'
#' @seealso \code{\link{ReadMarker}}, \code{\link{ReadPheno}}, and \code{\link{ReadMap}}
#'
#' @return
#' A list with the following components:
#' \describe{
#'\item{trait}{column name of the trait being used by 'AM'.}
#'\item{fformula}{Right hand size formula of the fixed effects part of the linear mixed model.}
#'\item{indxNA}{a vector containing the row indexes of those individuals, whose trait and fixed effects data contain
#' missing values and have been removed from the analysis.}
#' \item{Mrk}{a vector with the names of the snp in strongest and significant association with the trait.If no loci are found to be 
#' significant, then this component is \code{NA}.}
#' \item{Chr}{the chromosomes on which the identified snp lie.}
#' \item{Pos}{the map positions for the identified snp.}
#' \item{Indx}{the column indexes in the marker file of the identified snp.} 
#' \item{ncpu}{number of cpu used for the calculations.}
#' \item{availmemGb}{amount of RAM in gigabytes that has been set by the user.}
#' \item{quiet}{ boolean value of the parameter.}
#' \item{extBIC}{numeric vector with the extended BIC values for the loci  found to be in  significant association with the trait.}
#'}
#'
#' @examples
#'   \dontrun{ 
#'   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
#'   # However, the code can be run by the user. 
#'   #
#'
#'   #-------------------------
#'   #  Example  
#'   #------------------------
#'
#'   # read the map 
#'   #~~~~~~~~~~~~~~
#'   
#'   # File is a plain space separated text file with the first row 
#'   # the column headings
#'   complete.name <- system.file('extdata', 'map.txt', 
#'                                    package='Eagle')
#'   map_obj <- ReadMap(filename=complete.name) 
#'
#'   # read marker data
#'   #~~~~~~~~~~~~~~~~~~~~
#'   # Reading in a PLINK ped file 
#'   # and setting the available memory on the machine for the reading of the data to 8  gigabytes
#'   complete.name <- system.file('extdata', 'geno.ped', 
#'                                      package='Eagle')
#'   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
#'  
#'   # read phenotype data
#'   #~~~~~~~~~~~~~~~~~~~~~~~
#'
#'   # Read in a plain text file with data on a single trait and two covariates
#'   # The first row of the text file contains the column names y, cov1, and cov2. 
#'   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
#'   
#'   pheno_obj <- ReadPheno(filename=complete.name)
#'            
#'
#'  # Performing multiple-locus genome-wide association mapping with a model 
#'  #    with no fixed effects except for an intercept. 
#'  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'  
#'   res <- AM(trait = 'y',
#'                            fformula=c('cov1+cov2'),
#'                            map = map_obj,
#'                            pheno = pheno_obj,
#'                            geno = geno_obj, availmemGb=8)
#' }
#'
AM <- function(trait=NULL, 
               fformula  = NULL,
               availmemGb=8, 
               geno=NULL, 
               pheno=NULL, 
               map = NULL,
               ncpu=detectCores(),
               ngpu=0,
               quiet=TRUE,
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
 ## ngpu            number of gpu available for computation


 ## print tile
 .print_title()

 ngpu <- 0  ### NEED TO CHANGE THIS WHEN gpu implemented. 


 error.code <- check.inputs.mlam(ncpu=ncpu , availmemGb=availmemGb, colname.trait=trait, 
                     map=map, pheno=pheno, geno=geno )
 if(error.code){
   message("\n The Eagle function AM has terminated with errors.\n")
   return(NULL)
 }



 ## checking if map is present. If not, generate a fake map. 
 if(is.null(map)){
   if(!quiet ){
     message(" Map file has not been supplied. An artificial map is being created but this map is not used in the analysis. \n")
     message(" It is only used for the reporting of results. \n")
   }
   ## map has not been supplied. Create own map
   map <- data.frame(SNP=paste("M", 1:geno[["dim_of_ascii_M"]][2], sep=""), 
                     Chr=rep(1, geno[["dim_of_ascii_M"]][2]), 
                     Pos=1:geno[["dim_of_ascii_M"]][2])
  }

 ## check that the number of rows in the map file match the number of columns in the geno file
 if (geno[["dim_of_ascii_M"]][2] != nrow(map)){
   message(" Error: There is a differing number of loci read in by ReadMarker and ReadMap functions. \n")
   message("         The number of marker loci read in by ReadMarker() is ", geno[["dim_of_ascii_M"]][2], "\n")
   message("        The number of marker loci in  the marker map is  ", nrow(map), "\n") 
   message("\n AM has terminated with errors.\n")
   return(NULL)
 }


 ## check that the number of rows in the phenotype file match the number of rows in the geno file
 if (geno[["dim_of_ascii_M"]][1] != nrow(pheno)){
   message(" Error: There is a differing number  of rows read in by ReadMarker and ReadPheno functions. \n")
   message("         The number of rows read in by ReadMarker() is ", geno[["dim_of_ascii_M"]][1], "\n")
   message("        The number of rows  read in by ReadPheno is  ", nrow(map), "\n") 
   message("\n AM has terminated with errors.\n")
   return(NULL)
 }




 selected_loci <- NA
 new_selected_locus <- NA
 extBIC <- vector("numeric", 0)
 ## assign trait 
 trait <-  pheno[[trait]]




 ## Turn fformula  into class formula with some checks
if(!is.null(fformula)){
 if(fformula=="")  ## added for shiny
      fformula<-NULL
 }
 if(!is.null(fformula) ){
   if(length(grep("~", fformula))==0){
      if(length(fformula)==1){
          fformula <- as.formula(paste("~", fformula, sep="") )
      }  else {
          message(" fformula has ", length(fformula), " separate terms. It should be a single formula. \n") 
          message("\n AM has terminated with errors.\n")
          return(NULL)
      }
   } else {
    ## problem: formula should not contain ~
    message(" It looks like fformula contains a formula. \n")
    message(" If so, only the terms on the right hand side of the formula should be specified. \n")
    message(" Please remove the ~ from the formula. \n")
    message("\n AM has terminated with errors.\n")
    return(NULL)
  }  ## if length grep
 } ## end if(!is.null(fformula))


  ## check that terms in  formula are in pheno file
 if(!is.null(fformula)){
  res <- tryCatch(
     mat <- get_all_vars(formula=fformula, data=pheno) , 
     error = function(e) 
     {
         return(TRUE)
     }
  )
  if(!is.data.frame(res))
  {
   if(res){
      message(" fformula contains terms that are not column headings in the phenotype file. \n")
      message(" Check spelling and case of terms in fformula. \n")
      message("\n  AM has terminated with errors.\n ")
      return(NULL)
   }
  }
 }



 
 ## check for NA's in explanatory variables 
 ## If any, set individual's trait value to NA
 ## This means this individual will later be removed. 
 if(!is.null(fformula)){
    mat <- get_all_vars(formula=fformula, data=pheno)
    mat.of.NA  <- which(is.na(mat), arr.ind=TRUE)
  if(!is.null(dim(mat.of.NA)[1]) ){
     if(dim(mat.of.NA)[1]>0){
       trait[unique(mat.of.NA[,1])] <- NA
     }
  }
 }

 ## check for NA's in trait
 indxNA <- check.for.NA.in.trait(trait=trait)



 ## remove missing observations from trait
 if(length(indxNA)>0){
    trait <- trait[-indxNA]

    if(!quiet ){
     message(" The following rows are being removed from pheno due to missing data: \n")
     message(cat("             ", indxNA, "\n\n"))
    }

 }


## create a new M.ascii and Mt.ascii if length(indxNA) is non-zero 
## remove rows in M.ascii and columns in Mt.ascii of those individuals listed in indxNA 
if(length(indxNA)>0){
    res <- ReshapeM(fnameM=geno$asciifileM, fnameMt=geno$asciifileMt, indxNA=indxNA, dims=geno$dim_of_ascii_M)
    message(cat("new dimensions of reshaped M", res, "\n"))

     if(.Platform$OS.type == "unix") {
       geno$asciifileM <- paste(tempdir() , "/", "M.asciitmp", sep="")
     } else {
       geno$asciifileM <- paste( tempdir() , "\\", "M.asciitmp", sep="")
     }

     if(.Platform$OS.type == "unix") {
       geno$asciifileMt <- paste( tempdir() , "/", "Mt.asciitmp", sep="")
     } else {
       geno$asciifileMt <- paste( tempdir() , "\\", "Mt.asciitmp", sep="")
     }

    geno$dim_of_ascii_M <- res
}



 ## build design matrix currentX
 currentX <- .build_design_matrix(pheno=pheno, indxNA=indxNA, fformula=fformula, quiet=quiet )

 ## check currentX for solve(crossprod(X, X)) singularity
 chck <- tryCatch({ans <- solve(crossprod(currentX, currentX))},
           error = function(err){
            return(TRUE)
           })

  if(is.logical(chck)){
      if(chck){
        message(" There is a problem with the effects in fformula.\n")
        message(" These effects are causing computational instability. \n")
        message(" This can occur when there is a strong dependency between the effects.\n")
        message(" Try removing some of the effects in fformula. \n")
        message("\n  AM has terminated with errors.\n")
        return(NULL)
      }
  }


 ## Initialization
 continue <- TRUE
 itnum <- 1


 while(continue){
  message("\n\n Iteration" , itnum, ": Searching for most significant marker-trait association\n\n")
   ## based on selected_locus, form model matrix X
  currentX <- constructX(fnameM=geno[["asciifileM"]], currentX=currentX, loci_indx=new_selected_locus,
                          dim_of_ascii_M=geno[["dim_of_ascii_M"]],
                          map=map, availmemGb = availmemGb)  



    ## calculate Ve and Vg
    Args <- list(geno=geno,availmemGb=availmemGb,
                    ncpu=ncpu,selected_loci=selected_loci,
                    quiet=quiet)

    if(itnum==1){
        if(!quiet)
           message(" quiet=FALSE: calculating M %*% M^t. \n")
         MMt <- do.call(.calcMMt, Args)  


         if(!quiet)
             doquiet(dat=MMt, num_markers=5 , lab="M%*%M^t")
        invMMt <- chol2inv(chol(MMt))   ## doesn't use GPU
        gc()
    } 
    if(!quiet){
      message(" Calculating variance components for multiple-locus model. \n")
    }
    vc <- .calcVC(trait=trait, currentX=currentX,MMt=MMt, ngpu=ngpu) 
    gc()
    best_ve <- vc[["ve"]]
    best_vg <- vc[["vg"]]



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
     ARgs <- list(geno=geno,availmemGb=availmemGb, selected_loci=selected_loci,
                 MMt=MMt, invMMt=invMMt, best_ve=best_ve, best_vg=best_vg, currentX=currentX,
                 ncpu=ncpu, quiet=quiet, trait=trait, ngpu=ngpu)
      new_selected_locus <- do.call(.find_qtl, ARgs)  ## memory blowing up here !!!! 
     gc()
     selected_loci <- c(selected_loci, new_selected_locus)

   }  else {
     ## terminate while loop, 
     continue <- FALSE
   }  ## end if else


   itnum <- itnum + 1
   ## alternate stopping rule - if maxit has been exceeded.
    if(itnum > maxit){
         continue <- FALSE 
         .print_header()
         ## need to remove the last selected locus since we don't go on and calculate its H and extBIC 
         ## under this new model. 
         .print_final(selected_loci[-length(selected_loci)], map, extBIC)
         sigres <- .form_results(trait, selected_loci[-length(selected_loci)], map,  fformula, 
                     indxNA, ncpu, availmemGb, quiet,  extBIC )   
    }
 
  }  ## end while continue

if( itnum > maxit){
    .print_header()
    .print_final(selected_loci, map,  extBIC)
    sigres <- .form_results(trait, selected_loci, map,  fformula, 
                     indxNA, ncpu, availmemGb, quiet,  extBIC )   

} else {
    ## remove last selected_loci as for this locus, the extBIC went up
    if(length(selected_loci)>1){
        .print_header()
        .print_final(selected_loci[-length(selected_loci)], 
                     map, 
                     extBIC[-length(selected_loci)])
        sigres <- .form_results(trait, selected_loci[-length(selected_loci)], map,  fformula, 
                         indxNA, ncpu, availmemGb, quiet, 
                         extBIC[-length(selected_loci)] )   
    } else {
        .print_header()
        .print_final(selected_loci, map, extBIC)
        sigres <- .form_results(trait, selected_loci, map,  fformula, 
                         indxNA, ncpu, availmemGb, quiet, extBIC )   
   }  ## end inner  if(length(selected_locus)>1)
}  ## end if( itnum > maxit)

 
return( sigres )

} ## end AM










