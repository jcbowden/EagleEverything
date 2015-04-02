### To run multiple CPU
## salloc --ntasks-per-node=1  --ntasks-per-node=6 --mem=20gb --time=30:0 srun --pty bash

## To DO
## 
## read in M with rows as plants or  with rows as colums (as with  human data)
## how to get dimension of M when it is a binary file


##-------------------------------------
##  EMMA code 
##------------------------------------
emma.delta.ML.dLL.w.Z <-  function (logdelta, lambda, etas.1, xi.1, n, etas.2.sq) 
{
    t <- length(xi.1)
    delta <- exp(logdelta)
    etasq <- etas.1 * etas.1
    ldelta <- lambda + delta
    return(0.5 * (n * (sum(etasq/(ldelta * ldelta)) + etas.2.sq/(delta * 
        delta))/(sum(etasq/ldelta) + etas.2.sq/delta) - (sum(1/(xi.1 + 
        delta)) + (n - t)/delta)))
}

 emma.eigen.L.w.Z <- function (Z, K, complete = TRUE) 
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    eig <- eigen(K %*% crossprod(Z, Z), symmetric = FALSE, EISPACK = TRUE)
    return(list(values = eig$values, vectors = qr.Q(qr(Z %*% 
        eig$vectors), complete = TRUE)))
}


emma.eigen.R.w.Z <-  function (Z, K, X, complete = TRUE) 
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    n <- nrow(Z)
    t <- ncol(Z)
    q <- ncol(X)
    SZ <- Z - X %*% solve(crossprod(X, X)) %*% crossprod(X, Z)
    eig <- eigen(K %*% crossprod(Z, SZ), symmetric = FALSE, EISPACK = TRUE)
    if (is.complex(eig$values)) {
        eig$values <- Re(eig$values)
        eig$vectors <- Re(eig$vectors)
    }
    qr.X <- qr.Q(qr(X))
    return(list(values = eig$values[1:(t - q)], vectors = qr.Q(qr(cbind(SZ %*% 
        eig$vectors[, 1:(t - q)], qr.X)), complete = TRUE)[, 
        c(1:(t - q), (t + 1):n)]))
}

emma.delta.REML.dLL.w.Z <- function (logdelta, lambda, etas.1, n, t1, etas.2.sq) 
{
    t <- t1
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    etasq <- etas.1 * etas.1
    ldelta <- lambda + delta
    return(0.5 * (nq * (sum(etasq/(ldelta * ldelta)) + etas.2.sq/(delta * 
        delta))/(sum(etasq/ldelta) + etas.2.sq/delta) - (sum(1/ldelta) + 
        (n - t)/delta)))
}

emma.delta.REML.LL.w.Z <- function (logdelta, lambda, etas.1, n, t, etas.2.sq) 
{
    tq <- length(etas.1)
    nq <- n - t + tq
    delta <- exp(logdelta)
    return(0.5 * (nq * (log(nq/(2 * pi)) - 1 - log(sum(etas.1 * 
        etas.1/(lambda + delta)) + etas.2.sq/delta)) - (sum(log(lambda + 
        delta)) + (n - t) * logdelta)))
}


emma.MLE <- function (y, X, K, Z = NULL, ngrids = 100, llim = -10, ulim = 10, 
    esp = 1e-10, eig.L = NULL, eig.R = NULL) 
{
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    if (det(crossprod(X, X)) == 0) {
        warning("X is singular")
        return(list(ML = 0, delta = 0, ve = 0, vg = 0))
    }
    if (is.null(Z)) {
        if (is.null(eig.L)) {
            eig.L <- emma.eigen.L.wo.Z(K)
        }
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.wo.Z(K, X)
        }
        etas <- crossprod(eig.R$vectors, y)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
        Lambdas <- matrix(eig.R$values, n - q, m) + matrix(delta, 
            n - q, m, byrow = TRUE)
        Xis <- matrix(eig.L$values, n, m) + matrix(delta, n, 
            m, byrow = TRUE)
        Etasq <- matrix(etas * etas, n - q, m)
        LL <- 0.5 * (n * (log(n/(2 * pi)) - 1 - log(colSums(Etasq/Lambdas))) - 
            colSums(log(Xis)))
        dLL <- 0.5 * delta * (n * colSums(Etasq/(Lambdas * Lambdas))/colSums(Etasq/Lambdas) - 
            colSums(1/Xis))
        optlogdelta <- vector(length = 0)
        optLL <- vector(length = 0)
        if (dLL[1] < esp) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim, 
                eig.R$values, etas, eig.L$values))
        }
        if (dLL[m - 1] > 0 - esp) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim, 
                eig.R$values, etas, eig.L$values))
        }
        for (i in 1:(m - 1)) {
            if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] > 
                0) && (dLL[i + 1] < 0)) {
                r <- uniroot(emma.delta.ML.dLL.wo.Z, lower = logdelta[i], 
                  upper = logdelta[i + 1], lambda = eig.R$values, 
                  etas = etas, xi = eig.L$values)
                optlogdelta <- append(optlogdelta, r$root)
                optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root, 
                  eig.R$values, etas, eig.L$values))
            }
        }
    }
    else {
        if (is.null(eig.L)) {
            eig.L <- emma.eigen.L.w.Z(Z, K)
        }
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.w.Z(Z, K, X)
        }
        etas <- crossprod(eig.R$vectors, y)
        etas.1 <- etas[1:(t - q)]
        etas.2 <- etas[(t - q + 1):(n - q)]
        etas.2.sq <- sum(etas.2 * etas.2)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
        Lambdas <- matrix(eig.R$values, t - q, m) + matrix(delta, 
            t - q, m, byrow = TRUE)
        Xis <- matrix(eig.L$values, t, m) + matrix(delta, t, 
            m, byrow = TRUE)
        Etasq <- matrix(etas.1 * etas.1, t - q, m)
        dLL <- 0.5 * delta * (n * (colSums(Etasq/(Lambdas * Lambdas)) + 
            etas.2.sq/(delta * delta))/(colSums(Etasq/Lambdas) + 
            etas.2.sq/delta) - (colSums(1/Xis) + (n - t)/delta))
        optlogdelta <- vector(length = 0)
        optLL <- vector(length = 0)
        if (dLL[1] < esp) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim, 
                eig.R$values, etas.1, eig.L$values, n, etas.2.sq))
        }
        if (dLL[m - 1] > 0 - esp) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim, 
                eig.R$values, etas.1, eig.L$values, n, etas.2.sq))
        }
        for (i in 1:(m - 1)) {
            if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] > 
                0) && (dLL[i + 1] < 0)) {
                r <- uniroot(emma.delta.ML.dLL.w.Z, lower = logdelta[i], 
                  upper = logdelta[i + 1], lambda = eig.R$values, 
                  etas.1 = etas.1, xi.1 = eig.L$values, n = n, 
                  etas.2.sq = etas.2.sq)
                optlogdelta <- append(optlogdelta, r$root)
                optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root, 
                  eig.R$values, etas.1, eig.L$values, n, etas.2.sq))
            }
        }
    }
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    if (is.null(Z)) {
        maxva <- sum(etas * etas/(eig.R$values + maxdelta))/n
    }
    else {
        maxva <- (sum(etas.1 * etas.1/(eig.R$values + maxdelta)) + 
            etas.2.sq/maxdelta)/n
    }
    maxve <- maxva * maxdelta
    return(list(ML = maxLL, delta = maxdelta, ve = maxve, vg = maxva))
}


emma.delta.ML.LL.wo.Z <- function (logdelta, lambda, etas, xi) 
{
    n <- length(xi)
    delta <- exp(logdelta)
    return(0.5 * (n * (log(n/(2 * pi)) - 1 - log(sum((etas * 
        etas)/(lambda + delta)))) - sum(log(xi + delta))))
}





emma.eigen.L.wo.Z <- function (K) 
{
    eig <- eigen(K, symmetric = TRUE)
    return(list(values = eig$values, vectors = eig$vectors))
}

emma.eigen.R.wo.Z <-  function (K, X) 
{
    n <- nrow(X)
    q <- ncol(X)
    S <- diag(n) - X %*% solve(crossprod(X, X)) %*% t(X)
    eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric = TRUE)
    stopifnot(!is.complex(eig$values))
    return(list(values = eig$values[1:(n - q)] - 1, vectors = eig$vectors[, 
        1:(n - q)]))
}


emma.delta.ML.LL.w.Z <-  function (logdelta, lambda, etas.1, xi.1, n, etas.2.sq) 
{
    t <- length(xi.1)
    delta <- exp(logdelta)
    return(0.5 * (n * (log(n/(2 * pi)) - 1 - log(sum(etas.1 * 
        etas.1/(lambda + delta)) + etas.2.sq/delta)) - (sum(log(xi.1 + 
        delta)) + (n - t) * logdelta)))
}

 emma.delta.ML.LL.w.Z <- function (logdelta, lambda, etas.1, xi.1, n, etas.2.sq) 
{
    t <- length(xi.1)
    delta <- exp(logdelta)
    return(0.5 * (n * (log(n/(2 * pi)) - 1 - log(sum(etas.1 * 
        etas.1/(lambda + delta)) + etas.2.sq/delta)) - (sum(log(xi.1 + 
        delta)) + (n - t) * logdelta)))
}

emma.delta.ML.dLL.wo.Z <- function (logdelta, lambda, etas, xi) 
{
    n <- length(xi)
    delta <- exp(logdelta)
    etasq <- etas * etas
    ldelta <- lambda + delta
    return(0.5 * (n * sum(etasq/(ldelta * ldelta))/sum(etasq/ldelta) - 
        sum(1/(xi + delta))))
}



 emma.REMLE <-  function (y, X, K, Z = NULL, ngrids = 100, llim = -10, ulim = 10, 
    esp = 1e-10, eig.L = NULL, eig.R = NULL) 
{
    n <- length(y)
    t <- nrow(K)
    q <- ncol(X)
    stopifnot(ncol(K) == t)
    stopifnot(nrow(X) == n)
    if (det(crossprod(X, X)) == 0) {
        warning("X is singular")
        return(list(REML = 0, delta = 0, ve = 0, vg = 0))
    }
    if (is.null(Z)) {
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.wo.Z(K, X)
        }
        etas <- crossprod(eig.R$vectors, y)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
        Lambdas <- matrix(eig.R$values, n - q, m) + matrix(delta, 
            n - q, m, byrow = TRUE)
        Etasq <- matrix(etas * etas, n - q, m)
        LL <- 0.5 * ((n - q) * (log((n - q)/(2 * pi)) - 1 - log(colSums(Etasq/Lambdas))) - 
            colSums(log(Lambdas)))
        dLL <- 0.5 * delta * ((n - q) * colSums(Etasq/(Lambdas * 
            Lambdas))/colSums(Etasq/Lambdas) - colSums(1/Lambdas))
        optlogdelta <- vector(length = 0)
        optLL <- vector(length = 0)
        if (dLL[1] < esp) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim, 
                eig.R$values, etas))
        }
        if (dLL[m - 1] > 0 - esp) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim, 
                eig.R$values, etas))
        }
        for (i in 1:(m - 1)) {
            if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] > 
                0) && (dLL[i + 1] < 0)) {
                r <- uniroot(emma.delta.REML.dLL.wo.Z, lower = logdelta[i], 
                  upper = logdelta[i + 1], lambda = eig.R$values, 
                  etas = etas)
                optlogdelta <- append(optlogdelta, r$root)
                optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root, 
                  eig.R$values, etas))
            }
        }
    }
   else {
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.w.Z(Z, K, X)
        }
        etas <- crossprod(eig.R$vectors, y)
        etas.1 <- etas[1:(t - q)]
        etas.2 <- etas[(t - q + 1):(n - q)]
        etas.2.sq <- sum(etas.2 * etas.2)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
        Lambdas <- matrix(eig.R$values, t - q, m) + matrix(delta, 
            t - q, m, byrow = TRUE)
        Etasq <- matrix(etas.1 * etas.1, t - q, m)
        dLL <- 0.5 * delta * ((n - q) * (colSums(Etasq/(Lambdas * 
            Lambdas)) + etas.2.sq/(delta * delta))/(colSums(Etasq/Lambdas) + 
            etas.2.sq/delta) - (colSums(1/Lambdas) + (n - t)/delta))
        optlogdelta <- vector(length = 0)
        optLL <- vector(length = 0)
        if (dLL[1] < esp) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim, 
                eig.R$values, etas.1, n, t, etas.2.sq))
        }
        if (dLL[m - 1] > 0 - esp) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim, 
                eig.R$values, etas.1, n, t, etas.2.sq))
        }
        for (i in 1:(m - 1)) {
            if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] > 
                0) && (dLL[i + 1] < 0)) {
                r <- uniroot(emma.delta.REML.dLL.w.Z, lower = logdelta[i], 
                  upper = logdelta[i + 1], lambda = eig.R$values, 
                  etas.1 = etas.1, n = n, t1 = t, etas.2.sq = etas.2.sq)
                optlogdelta <- append(optlogdelta, r$root)
                optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root, 
                  eig.R$values, etas.1, n, t, etas.2.sq))
            }
        }
    }
    maxdelta <- exp(optlogdelta[which.max(optLL)])
    maxLL <- max(optLL)
    if (is.null(Z)) {
        maxva <- sum(etas * etas/(eig.R$values + maxdelta))/(n - 
            q)
    }
    else {
        maxva <- (sum(etas.1 * etas.1/(eig.R$values + maxdelta)) + 
            etas.2.sq/maxdelta)/(n - q)
    }
    maxve <- maxva * maxdelta
    return(list(REML = maxLL, delta = maxdelta, ve = maxve, vg = maxva))
}


emma.delta.REML.dLL.wo.Z <-  function (logdelta, lambda, etas) 
{
    nq <- length(etas)
    delta <- exp(logdelta)
    etasq <- etas * etas
    ldelta <- lambda + delta
    return(0.5 * (nq * sum(etasq/(ldelta * ldelta))/sum(etasq/ldelta) - 
        sum(1/ldelta)))
}

emma.delta.REML.LL.wo.Z <-  function (logdelta, lambda, etas) 
{
    nq <- length(etas)
    delta <- exp(logdelta)
    return(0.5 * (nq * (log(nq/(2 * pi)) - 1 - log(sum(etas * 
        etas/(lambda + delta)))) - sum(log(lambda + delta))))
}






#### To run multple GPU's
### > export OMP_NUM_THREADS=$PBS_NUM_PPN
## > module load cuda/6.0 R
## > module load R/3.0.0
## > LD_PRELOAD=libnvblas.so R  
## monitoring gpu usage
##    nvidia-smi -l 3


##library('Rcpp')
##library('RcppEigen')
##library('matrixcalc')
##library('Matrix')

##---------------------------
## Rcpp Function Declarion
##---------------------------

## This builds a dll for the function
## sourceCpp("RcppFunctions.cpp", rebuild=TRUE, verbose=TRUE)


##-------------------------------
## R Function Declaration 
##-------------------------------
#' @title Check correctness of marker genotypes.
#' 
#' @description
#' \code{check.genofile} performs various checks on the marker genotype file.
#' @param fnameIN  character vector containing the name of the marker genotype file
#' @param dirPath  character vector contain the directory path to where the marker genotype file is located. 
#' @param columnwise a logical value. When \code{TRUE},  each row contains an individuals genotypes and each 
#'           column contains a marker locus\' genotypes. When \code{FALSE}, each row contains a marker locus\'s
#'           genotypes and each  column contains an individual\'s genotypes. 
#'                    When \code{FALSE}, a row contains the data from a marker locus. 
#' @param AA       integer value corresponding to the AA genotype in the marker genotype file. This must be specified. 
#' @param AB       integer value corresponding to the AB genotype in the marker genotype file. This can be left unspecified 
#'                  if there are no hets.
#' @param BB       integer value corresponding to the BB genotype in the marker genotype file. This must be specified. 
#' @param check_num_geno_in_row a logical value. When \code{TRUE}, the number of genotypes in each row is returned.
#' @param check_genotypes a logical value. When \code{TRUE}, it checks that the marker genoyptes are either 0, 1, or 2. 
#' @param csv      is a logical value. When \code{TRUE}, a comma seperated version (csv) file is assumed. 
#'
#' @details
#' The function \code{check.genofile} checks the following that:
#' \itemize{
#'   \item{ the file exists}
#'   \item{ there is the same number of genotypes per row}
#'   \item{  the file only contains numeric values as specified by the parameters AA, AB, and BB.}a
#'   \item{  the file only contains two (if \code{AB} has not been set) or three different numeric values.}
#' }
#'
#'  Missing values are not permitted and will cause an error. 
#'
#' The marker genotype file must have the marker genotypes coded as numeric values. Any numeric values are valid. However, 
#' the parameters \code{AA} and \code{BB} must be specified. If the parameter \code{AB} is not set, an inbred population 
#' is assumed and any  heterozygous genotypes found in the marker file will cause an error. 
#' @return  a list is returned with elements \code{file_exists} and \code{num_genotypes_per_row}. \code{num_genotypes_per_row}
#'          will be \code{NULL} unless the \code{check_num_geno_in_row} parameter has been set to \code{TRUE}.
#' @seealso \code{\link{read.genotypes}}
check.genofile <- function(fnameIN=NULL, dirPath=getwd(),
                           columnwise= TRUE, 
                           AA=NULL,
                           AB=NULL,
                           BB=NULL,
                           check_num_geno_in_row=FALSE,
                           check_genotypes=FALSE, 
                           csv = FALSE)
{
 ## function to perform a series of tests on the genotype file
 ## Args
 ##      fnameIN     file name of ASCII genotype file (no header file handling yet)
 ##      dirPath     path to ASCII genotype file
 ##      check_num_geno_in_row  if TRUE, then count.fields is run to check that the 
 ##                             number of fields per row is the same. This can take a 
 ##                             long time if the file is large. 
 ##      check_genotypes        if TRUE, each row is checked for genotypes that are not 0,1,2

  res <- list(file_exists=TRUE, num_genotypes_per_row=NULL) 


  if(.Platform$OS.type == "unix") {
    dir_path  <- paste(dirPath, "/",  sep="")
  } else {
    dir_path    <- paste(dirPath, "\\",  sep="")
  }

  file_geno <- paste(dir_path, fnameIN, sep="")



  ## does the file exist
  if(!file.exists(file_geno))
  {
    res[["file_exists"]] <- FALSE
    cat(" File ", file_geno, " does not exist. \n")
    stop(" Please modify path and/or name of genotype file. \n")
  }


  ## has genotype  file name been provided
  if(is.null(fnameIN))
      stop(" No name of the genotype file has been supplied.")


  ## does the file contain a constant number of records. 
  if(check_num_geno_in_row){
     res[["num_genotypes_per_row"]] <- count.fields(file=file_geno)
     if(length(unique( res[["num_genotypes_per_row"]] ))>1){
       cat(" File ", file_geno, " has  differing numbers of genotypes per line ... \n")
       cat(" The differing number of gentoypes per line are ", unique(res[["num_genotypes_per_row"]]), "\n")
       stop(" Please modify the genotype file to have the same number of genotypes per row. \n")
     } else {
       cat(" File ", file_geno, " has ", unique( res[["num_genotypes_per_row"]] ), " genotypes per line...\n\n")
     } ## end if else
  }  ## end if checkrownum


  ## Has AA, AB, BB been assigned numeric values
  if(is.null(AA) |  is.null(BB))
  {
     stop(" AA and BB must be assigned a numeric value \n")
  }

  if(!is.numeric(AA) | !is.numeric(BB))
     stop(" AA and/or BB must be an integer value corresponding to the AA and/or BB genotype, respectively,  in the marker genotype file. ")
  if(!is.null(AB))
     if(!is.numeric(AB))
       stop(" AB must be an integer value corresponding to the AB genotype in the marker genotype file. ")



  ## does the file contain 0,1,2 genotypes only (implemented in Rcpp for speed)
  if(check_genotypes){
     if(is.null(AB))
         AB <- 18923282  ## no hets so setting it to something weird
     checkGenotypes(file_geno, AA, AB, BB, csv)
}
  return(res)
}






calculateMMt <- function(geno=NULL, workingmemGb, numcores, selected_loci=NA, dim_of_bin_M=NULL)
{
 ## R interface to Rcpp code to calculate M %*% t(M)
 ## Args
 ##      geno        absolute path + file name of binary packed M file
 ##      workingmemGb    amount of memory in Gbytes available for creation of MMt
 ##      numcores    number of cores for matrix operations
 ##      selectedloci an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                   have been selected to act as fixed QTL effects in the model. 
 ##      dim_of_bin_M    numeric vector with the row, column numbers of M. 
  #------------------------------------------
  # bin file about to be overwritten
  #------------------------------------------

  if(!file.exists(geno)){
    cat(" Error: The binary packed file ", geno, " cannot be found.\n")
    stop(" calculateMMt has terminated with errors.") 
   }


  MMt <- calculateMMt_rcpp( f_name_bin=geno, selected_loci = selected_loci,
                               max_memory_in_Gbytes=workingmemGb, num_cores=numcores, 
                               dims= dim_of_bin_M)
  return(MMt)

}  ## end function









calculateMMt_sqrt_and_sqrtinv <- function(MMt=NULL, checkres=TRUE, verbose=FALSE)
{
  ## R function for calculating the square root of M * M^t
  ## and the inverse of the square root of MMt
  ## where M * M^t has already been created. 
  ## Using SVD for calculation
  ## Args
  ##  MMt   a matrix of M * M^t
  ##  checkres  when true, the accuracy of the inversion is checked. 

  ## testing that MMt is postive definite
  if(!is.positive.definite(MMt))
    stop(" The MMt matrix is not postive definite. This can occur if you have duplicate rows in  the genotype file.")
   
   ## calculate square root of MMt 
   cat(" WARNING: this may take some time if the genotype file is large and/or computations \n")
   cat("  are not being distributed across multiple threads.\n\n")
   cat(" Beginning SVD calculation ... \n")
   svdM <- svd(MMt)
   Sigma <- diag(sqrt(svdM[["d"]]))
   sqrt_MMt <- svdM[["u"]] %*% Sigma %*% t(svdM[["v"]])
   rm(Sigma)
   gc()

   ## calculate inverse square root of MMT
   inverse_Sigma <- diag(1/sqrt(svdM[["d"]]))
   inverse_sqrt_MMt <- svdM[["u"]] %*% inverse_Sigma %*% t(svdM[["v"]])
   rm(inverse_Sigma)
   gc()
   if(checkres){
       a <- (sqrt_MMt %*% inverse_sqrt_MMt)
       if(trunc(sum(diag(a))) != nrow(MMt))
       {
         cat(" \n\n\nWARNING: these results may be unstable.\n")
         cat(" The sum of the diagonal elements of the square root of MMt and its inverse is ", sum(diag(a)), " where \n")
         cat("  it should have been ", nrow(MMt), "\n")
         cat("  This can occur if the genotype file contains near identical rows.  Please check.\n\n")
      
         if(verbose)
         {
           cat(" --- diagonal  values follow --- \n")
           cat( diag(a))
           cat(" Diagonal elements different from 1 indicates that the genotype file contains individuals with near identical marker genotypes. \n")
         }


       } 
   }   ## end if(checkres)
   res <- list(sqrt_MMt=sqrt_MMt, inverse_sqrt_MMt=inverse_sqrt_MMt)
  


} ## end function



calculateH <- function(MMt=NULL, varE=NULL, varG=NULL)
{
  ## R function for calculating the H variance matrix 
  ## which is
  ##  H = \sigma^2_E I  + \sigma^2_G  MMt
  ## Args:
  ##     MMt  - matrix object for M %*% M^T
  ##     varE  -  numeric value for the residual variance
  ##     varG  -  numeric value for the polygenic variance (\sigma^2_g)
  ##
  ## H matrix is returned. 

  if(!is.numeric(varE))
    stop(" The varE (residual variance) must be numeric.")

  if(varE < 0)
    stop(" VarE cannot be negative.")

  if(!is.numeric(varG))
    stop(" The varG (genotypic variance) must be numeric.")

  if(varG < 0)
    stop(" VarG cannot be negative.")


  if(is.null(MMt))
    stop("MMt cannot be null.")


  return( varE * diag(nrow(MMt)) + varG * MMt)


}


calculateP  <- function(H=NULL, X=NULL)
{
  ## R function to calculate P matrix
  ## Args:
  ##       H is the variance matrix
  ##       X is the design matrix supplied by the user
  ## Returns:
  ##   matrix object P

  if(is.null(H))
    stop(" H must be specified.")
  if(is.null(X))
    stop(" A design matrix has not be specified. ")

   if(nrow(H) != nrow(X))
      stop(" The number of rows in H and X are not the same.")


   Hinv <- solve(H)
   P <- Hinv - Hinv %*% X %*% solve( t(X) %*% Hinv %*% X )  %*% t(X) %*% Hinv


  return(P)

}


calculate_reduced_a <- function(varG=NULL, P=NULL, MMtsqrt=NULL, y=NULL)
{

  if( !(nrow(P) ==  length(y))){
    cat(" Error:  there is a problem with the  dimensions of  P, and/or the vector y.")
    cat("         They should  be of the dimension (n x n), and a vector of length n.")
    cat(" The dimensions are: \n")
    cat(" dim(P)      = ", dim(P), "\n")
    cat(" length(y)   = ", length(y), "\n")
    stop()

  }

 if(is.null(varG))
   stop(" VarG must be specified.")

  if(is.null(P))
   stop(" P must be specified")


  if(is.null(y))
   stop(" y must be specified")


    a <- varG * MMtsqrt %*% P %*% y

return(a)

}



mistake_calculate_reduced_a <- function(varG=NULL, bin_path=getwd(), P=NULL, y=NULL, workingmemGb=8, dim_of_bin_M=NULL, 
                                 selected_loci=NA)
{
 ## Rcpp function to calculate the BLUP (a) values under a dimension reduced model
 ## Args:
 ##    varG is a scalar value
 ##    bin_path  is bin director to Mt which is a p x n matrix
 ##    P   is a n x n matrix
 ##    y   is a n x 1 vector
 ##
 ## a* = sigma^2_a * t(M) * P * y
 ## Returns:
 ##   a numeric vector of dimension reduced a values 

 if(is.null(varG))
   stop(" VarG must be specified.")

  if(is.null(P))
   stop(" P must be specified")

 
  if(is.null(y))
   stop(" y must be specified")

 
  if( !(nrow(P) ==  length(y))){
    cat(" Error:  there is a problem with the  dimensions of  P, and/or the vector y.")
    cat("         They should  be of the dimension (n x n), and a vector of length n.")
    cat(" The dimensions are: \n")
    cat(" dim(P)      = ", dim(P), "\n")
    cat(" length(y)   = ", length(y), "\n")
    stop()  

  }

if(.Platform$OS.type == "unix") {
    bin_path  <- paste(bin_path, "/",  sep="")
} else {
   bin_path    <- paste(bin_path, "\\",  sep="")
}


  ycolmat <- matrix(data=y, ncol=1)  ## makes it easier when dealing with this in Rcpp
  fnamebin <- paste(bin_path, "Mt.bin", sep="")
  ar <- calculate_reduced_a_rcpp(f_name_bin = fnamebin, varG=varG, P=P, y=ycolmat, max_memory_in_Gbytes=workingmemGb, 
                                 dims=dim_of_bin_M , selected_loci = selected_loci )





## t(t(y)) is a trick to get y as a row matrix 
##return( varG * invMMt %*% P %*% t(t(y)))
return(ar)

}





calculate_a_and_vara <- function(bin_path=getwd(), maxmemGb=8, dims=NULL,
                         selectedloci = NA,
                         invMMtsqrt=NULL, transformed_a=NULL, transformed_vara=NULL)
{
 ## an Rcpp function to take dimension reduced a (BLUP) values 
 ## and transform them into the original a (BLUP) values and their variances 
 ## Args:
 ##   bin_path         path to the location of the binary file containing the 
 ##                              transposed M matrix
 ##   maxmemGb         maximum available memory (in Gigabytes) that are available for use
 ##   dims             a 2 element numeric vector with the number of rows,columns in M 
 ##   invMMtsqrt       a matrix object of the form (M %*% M^T)^{-0.5}
 ##   transformed_a    a numeric vector of the dimension reduced BLUP or a values
 ##   transformed_vara a numeric matrix of dimension dims(1) x dims(1) for the dimension reduced BLUPs (or a) values. 
 ##   selectedloci     an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                    have been selected to act as fixed QTL effects in the model. 



if(.Platform$OS.type == "unix") {
    bin_path  <- paste(bin_path, "/",  sep="")
} else {
   bin_path    <- paste(bin_path, "\\",  sep="")
}



  file_bin <- paste(bin_path, "Mt.bin",sep="")

  if(!file.exists(file_bin)){
      cat("\n\n  ERROR: ", file_bin, " does not exist and it should have been created. \n\n")
      stop()
  }

  dimsMt <- c(dims[2], dims[1]) 
  calculate_a_and_vara_rcpp(f_name_bin=file_bin,max_memory_in_Gbytes=maxmemGb, dims=dimsMt, 
                    selected_loci = selectedloci,
                    inv_MMt_sqrt=invMMtsqrt,  
                    dim_reduced_vara = transformed_vara,
                    a = transformed_a)

}


calculate_reduced_vara <- function(X=NULL, varE=NULL, varG=NULL, invMMt=NULL, MMtsqrt=NULL)
{
## Using var(\hat(a)) = simgaG - Cjj  where Cjj is the component from C^-1 (henerdsons 
##   mixed model equations coefficent matrix.   See Verbyla et al. TAG 2007.

##  Mixed model equations for the linear mixed model
##
##  X^T %*% R^-1  X                  X^T %*% R^-1 %*% Ze
##
##
##  Ze^t %*% R^-1 %*% X            Ze^t %*% R^-1 %*% Ze   +  G^-1
##
##  Ze = MMt^0.5
##  R  = (varE * I)^-1
##  G  = (varG * I)^-1
## 

  ## first principals
  Ze <- MMtsqrt
  R1  <- solve( varE * diag(nrow(invMMt)))
  G1  <- solve( varG * diag(nrow(invMMt)))
  A <- t(X) %*% R1 %*% X
  B <- t(X) %*% R1 %*% Ze
  C <- t(Ze) %*% R1 %*% X
  D <- t(Ze) %*% R1 %*% Ze + G1

  D1 <- solve(D)


#  C1 <- cbind(A,B)
#  C2 <- cbind(C,D)
#  CC <- rbind(C1, C2)
#  invCC <- solve(CC)

  ## ?? not sure about this ....
  vars <- varG * diag(nrow(D1))  - ( D1 + D1 %*% C %*% solve(A - B %*% D1 %*% C) %*% B %*% D1 )

    return(vars )

}



check.inputs <- function(numcores=NULL, workingmemGb=NULL, path=NULL, 
                         bin_path=NULL, file_genotype=NULL, file_phenotype=NULL, alpha=NULL){


if(!is.null(alpha))
{
  if(!is.numeric(alpha))
    stop(" alpha parameter is of wrong class. It should be a numeric.")

  if(alpha < 0 | alpha > 0.5)
    stop(" alpha value must be in the range 0-0.5.")

}


if(!is.null(numcores)){
 if(!is.numeric(numcores))
   stop(" numcores is not a numeric.")

 if(numcores < 1)
    stop(" numcores cannot be zero or a negative number.")
}

if(!is.null(workingmemGb))
{
 if(!is.numeric(workingmemGb))
   stop(" workingmemGb is not a numeric.")

 if(workingmemGb <= 0)
    stop(" workingmemGb cannot be zero or a a negative number.")
}

if(!is.null(path))
  if(!file.exists(path))
    stop(" path directory does not exist.")


if(!is.null(bin_path))
  if(!file.exists(bin_path))
      stop(" bin_path directory does not exist.")

if(!is.null(file_genotype))
{
  if(.Platform$OS.type == "unix") {
    genofile <- paste(path, "/", file_genotype, sep="")
  } else {
   genofile <- paste(path, "\\", file_genotype, sep="")
  }
  if(!file.exists(genofile)){
    msg <- paste(" Cannot find genotype file ", genofile, sep="")
    stop(msg)
  }
}


if(!is.null(file_phenotype))
{ 
  if(.Platform$OS.type == "unix") {
    phenofile <- paste(path, "/", file_phenotype, sep="")
  } else {
   phenofile <- paste(path, "\\", file_phenotype, sep="")
  }

  if(!file.exists(phenofile)){
    msg <- paste(" Cannot find phenotype file ", phenofile, sep="")
    stop(msg)
  }
}




}

#' @title Read phenotype file
#' @description Read in the phenotypic data
#' @param path  a character vector containing the absolute path for where the phenotype file is located.
#' @param file_phenotype a character vector containgin the name of the phenotype file.
#' @param header a logical value. When \code{TRUE}, the first row of the file contains the names of the columns. 
#' @param csv   a logical value. When \code{TRUE}, a csv file format is assumed. When \code{FALSE}, a space separated format is assumed. 
#' @details
#' Here, \code{read.phenotypes} reads in phenotypic data into the package. A space  separated
#' ASCII file with column headings is assumed. A csv file can also be read if \code{csv} is set to \code{TRUE}, 
#' Missing values are allowed and should be coded \code{NA}. However, if covariates contain 
#' missing values, and these covariates are used by \code{\link{multiple_locus_am}}, then  
#' \code{\link{multiple_locus_am}} will return an error.  \code{NA} values are allowed in the 
#' trait values. See \code{\link{multiple_locus_am}} for details. 
#' @seealso \code{\link{read.genotypes}}
#' @return 
#' a data frame is returned of the phenotypic data. If \code{header} is true, the 
#' names of the columns will be as specified by the first row of the phenotypic file. If \code{header} is \code{FALSE}, 
#' generic names are supplied by R in the form of V1, V2, etc.   
#'
#' @examples
#' # Read in  example phenotypic data from ./inst/extdata/
#' 
#' # find the full location of the phenotypic data 
#' complete.name <- system.file("extdata", "phenoexample.csv", package="WMAM")
#'   
#' # read in phenotypic data which is in csv format
#' phenodata <- read.phenotypes(path=dirname(complete.name),  
#'                              file_phenotype=basename(complete.name), 
#'                              csv=TRUE) 
#'                                
#'  ## print a couple of lines of the data file
#'  head(phenodata)
read.phenotypes <- function(path=getwd() , file_phenotype = NULL, header=TRUE, csv=FALSE){

  check.inputs(path=path, file_phenotype=file_phenotype)

  if(.Platform$OS.type == "unix") {
    phenofile <- paste(path, "/", file_phenotype, sep="")
  } else {
   phenofile <- paste(path, "\\", file_phenotype, sep="")
  }

  if(csv) sep=","
  phenos <- read.table(phenofile, header=header, sep=sep)
cat("\n\n Reading Phenotype File \n\n")
cat(" Loading file ....... \n\n")
cat("                    SUMMARY OF PHENOTYPE FILE  \n")
cat(" file location(path):         ",  dirname(phenofile), "\n")
cat(" file name:                   ",  basename(phenofile), "\n")
cat(" number of rows:              ", nrow(phenos), "\n")
cat(" number of columns:           ", ncol(phenos), "\n")
cat("                    Column classes  \n")
for(ii in 1:ncol(phenos))
  cat(c( sprintf("%20s   %15s", names(phenos)[ii], class(phenos[[ii]]) ), "\n"))


cat(" Warning!  If the column classes are incorrect, these will need to be changed by the user.")

  return(phenos)


}







#' @title Read map file
#' @description Read in the marker map  data
#' @param path  a character vector containing the absolute path for where the map file is located.
#' @param file_map a character vector containing the name of the map file.
#' @param csv   a logical value. When \code{TRUE}, a csv file format is assumed. When \code{FALSE}, a space separated format is assumed. 
#' @details
#' Here, \code{read.map} reads in map data into the package. A space  separated
#' ASCII file with column headings is assumed. A csv file can also be read if \code{csv} is set to \code{TRUE}, 
#' Missing values are not allowed. 
#' 
#' The order of the columns must be 
#'   column 1:  marker name
#'   column 2:  chromosome number
#'   column 3:  map position 
#'
#' The order of the markers in this file must match the column order of the genotype file because 
#' the genotype file does not contain marker names. 
#' @seealso \code{\link{read.genotypes}}
#' @return 
#' a data frame is returned of the map data. 
#'
#' @examples
#' # Read in  example map data from ./inst/extdata/
#' 
#' # find the full location of the map data 
#' complete.name <- system.file("extdata", "mapexample.txt", package="WMAM")
#'   
#' # read in map data 
#' mapdata <- read.map(path=dirname(complete.name),  
#'                              file_map=basename(complete.name)) 
#'                                
read.map  <- function(path=getwd(), file_map = NULL, csv=FALSE)
{
 check.inputs(path=path, file_phenotype=file_map)

  if(.Platform$OS.type == "unix") {
    mapfile <- paste(path, "/", file_map, sep="")
  } else {
   mapfile <- paste(path, "\\", file_map, sep="")
  }
  sep=""
  if(csv) sep=","
  map <- read.table(mapfile, header=TRUE, sep=sep)
cat("\n\n Reading Map File \n\n")
cat(" Loading file ....... \n\n")
cat("                    SUMMARY OF MAP FILE  \n")
cat(" file location(path):         ",  dirname(mapfile), "\n")
cat(" file name:                   ",  basename(mapfile), "\n")
cat(" number of markers:           ", nrow(map), "\n")
cat(" number of columns:           ", ncol(map), "\n")
cat(" number of chromosomes:       ", length(unique(map[[2]])), "\n")
cat(" first 10 markers of the map file ... \n")
print(head(map, n=10))
cat("\n\n")

return(map)

}



create.bin  <- function(file_genotype, bin_path, columnwise, AA, AB, BB, workingmemGb, dim_of_bin_M, csv){
 ## an Rcpp function to create the packed binary file of the genotype data M and Mt
 ## from genotype data that may be saved row or column wise in terms of the markers.
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## bin_path         path name for where binary files are to be stored. 
 ## columnwise       logical value. If TRUE, then the cols of the  genotype data are the markers. 
 ## AA, AB, BB       numeric codes for associated genotypes in marker genotype file
 ## workingmemGb     available memory for converstion to packed binary
 ## dim_of_bin_M             row, column dimensions of M.  

 if(columnwise){
    ## Here, M.bin can be created without transposing
    binMfile <- paste(bin_path, "M.bin", sep="") ## file name for binary packed Mt file.

    createM_rcpp(f_name = file_genotype, f_name_bin = binMfile, AA = AA, AB = AB, BB = BB,
               max_memory_in_Gbytes=workingmemGb,  dims = dim_of_bin_M , csv = csv)

    ## Mt.bin is created by transposing ASCII data
    binMtfile <- paste(bin_path, "Mt.bin", sep="") ## file name for binary packed Mt file.

    createMt_rcpp(f_name = file_genotype, f_name_bin = binMtfile,  AA = AA, AB = AB, BB = BB,
                  max_memory_in_Gbytes=workingmemGb,  dims = dim_of_bin_M, csv=csv )

 } else {
   ## Here, Mt.bin can be created without transposing
   binMtfile <- paste(bin_path, "Mt.bin", sep="") ## file name for binary packed Mt file.

   createM_rcpp(f_name = file_genotype, f_name_bin = binMtfile, AA = AA, AB = AB, BB = BB,
               max_memory_in_Gbytes=workingmemGb,  dims = dim_of_bin_M , csv = csv)

    ## Here, M.bin is created by transposing data
    binMfile <- paste(bin_path, "M.bin", sep="") ## file name for binary packed Mt file.

    createMt_rcpp(f_name = file_genotype, f_name_bin = binMfile,  AA = AA, AB = AB, BB = BB,
                  max_memory_in_Gbytes=workingmemGb,  dims = dim_of_bin_M, csv=csv )
 }

 return(NULL)

}

create.bin.Mt <- function(file_genotype, bin_path, columnwise, AA, AB, BB, workingmemGb, dim_of_bin_M, csv){
 ## an Rcpp function to create the packed binary file of the genotype data Mt
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## bin_path         path name for where binary files are to be stored. 
 ## AA, AB, BB       numeric codes for associated genotypes in marker genotype file
 ## workingmemGb     available memory for converstion to packed binary
 ## dim_of_bin_M             row, column dimensions of M.  

 binMtfile <- paste(bin_path, "Mt.bin", sep="") ## file name for binary packed Mt file.

 createMt_rcpp(f_name = file_genotype, f_name_bin = binMtfile,  AA = AA, AB = AB, BB = BB, 
               max_memory_in_Gbytes=workingmemGb,  dims = dim_of_bin_M, csv=csv )

 return(NULL)

}



create.bin.M <- function(file_genotype, bin_path, columnwise, AA, AB, BB, workingmemGb, dim_of_bin_M, csv){
 ## an Rcpp function to create the packed binary file of the genotype data M
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## bin_path         path name for where binary files are to be stored. 
 ## AA, AB, BB       numeric codes for associated genotypes in marker genotype file
 ## workingmemGb     available memory for converstion to packed binary
 ## dim_of_bin_M             row, column dimensions of M.  

 binMfile <- paste(bin_path, "M.bin", sep="") ## file name for binary packed Mt file.

 createM_rcpp(f_name = file_genotype, f_name_bin = binMfile, AA = AA, AB = AB, BB = BB, 
               max_memory_in_Gbytes=workingmemGb,  dims = dim_of_bin_M , csv = csv)

 return(NULL)


}



#' @title Read  marker genotype file.
#' 
#' @description
#' \code{read.genotypes} reads in the genotypes from the marker genoytpe file.
#' @param path  character vector containing the directory path to where the marker genotype file is located.
#' @param bin_path  character vector containing the directory path to where the binary converted marker genotype files are to be located. 
#' @param columnwise a logical value. When \code{TRUE},  each row contains an individuals genotypes and each 
#'           column contains a marker locus\' genotypes. When \code{FALSE}, each row contains a marker locus\'s
#'           genotypes and each  column contains an individual\'s genotypes. 
#'                    When \code{FALSE}, a row contains the data from a marker locus. 

#' @param AA       integer value corresponding to the AA genotype in the marker genotype file. This must be specified. 
#' @param AB       integer value corresponding to the AB genotype in the marker genotype file. This can be left unspecified  if there are no hets. 
#' @param BB       integer value corresponding to the BB genotype in the marker genotype file. This must be specified. 
#' @param file_genotype character vector containing the name of the marker genotype file.
#' @param  check  a logical value. If \code{TRUE}, then \code{\link{check.genofile}} is run to check the marker genotype file
#'         for errors. 
#' @param workingmemGb a numeric value. It specifies the amount of memory (in Gigabytes) available for reading in the marker 
#'                     genotype file. 
#' @param csv   a logical value. When \code{TRUE}, a csv file format is assumed. When \code{FALSE}, a space separated format is assumed. 
#'
#' @details
#'
#' \code{read.genotypes} reads the marker genoytpes from file and converts it into two binary files:
#' a binary M file and a binary Mt file. The binary M file is where the rows are the individuals and 
#' the columns are the marker loci.  The binary Mt file is where the rows are the marker loci and the columns
#' are the individuals.  This is true regardless of the value of \code{columnwise}. 
#'
#' \code{columnwise} allows genotype data to be read where the data may be organized as 
#' each column contains the genotypes for a marker locus (\code(columnwise)=\code{TRUE}) or 
#' each row contains the genotypes for a marker locus (\code{columnwise)=\code{FALSE}).
#'
#' The two binary files created by \code{read.genotypes} are called M.bin and Mt.bin and they are created in 
#' directory \code{bin_path}. 
#'
#' It is assumed that the marker genotype file is a space separated ASCII file but csv files can also 
#' be read by setting \code{csv} to \code{TRUE}.
#'
#' Missing values are not allowed. If present, this will cause and error. Such values should be replaced by their 
#' imputed values. 
#'
#' The \code{path} and \code{bin_path} parameters, if not set, will default to the current directory. 
#'
#' The parameters \code{AA} and \code{BB} must be specified to their associated numeric values in the marker 
#' genotype file.  If there are no heterozygous genotypes (i.e., it is an inbred population), 
#' then \code{AB} need not be specified.
#'
#' The \code{workingmemGb} parameter should be set to the largest amout of RAM for the machine, when practical.  If the marker 
#' genotype file is too large for the available memory, then it converted, chunk-wise,  into a single binary packed file.  
#' The consequences of reading the marker genotype file chunk- or block-wise is a longer processing time. 
#'
#' @return  a list is returned with elements \code{binfileM} , \code{binfileMt}, and \code{dim_of_bin_M}. 
#' which is the absolute path and file name of the binrary file for the marker genotype data, 
#' the absolute path and file name of the binary file for the transpose of the marker genotype data, 
#' and the number of rows and columns in the binary M genotype file, respectively.
#' 
#' @examples
#'   # find the full location of the genotype data. Data contained in ./inst/extdata/. 
#'   complete.name <- system.file("extdata", "genoexample.txt", package="WMAM")
#'
#'   
#'   # read in the ASCII marker genotype data where 0 values are being treated as genotype AA 
#'   # and 1 values are being treated as genoytpe BB. There are no heterozygotes so AB is not specified. 
#'   # 3 Gbytes of memory has been specified. The file is space separated with the rows the plants
#'   # and the columns the snp loci.
#'   geno.list <- read.genotypes(path=dirname(complete.name), AA=0, BB=1, 
#'                  file_genotype=basename(complete.name),  workingmemGb=2) 
#'    
#'
#'   # geno.list is a list with elements binfileM, binfileMt, and dim_of_bin_M
#'   # which corresponds to the name of the binary packed file for the marker genotype data
#'   # the name of the binary packed file for the transpose of the marker genotype data,
#'   # and a vector containing the number of rows and columns in the marker gentoype file. 
#'   (geno.list)
#'
#' @seealso \code{\link{check.genofile}}
read.genotypes <- function(path=getwd(), bin_path=getwd(), columnwise=TRUE, 
                           AA=NULL, AB=NULL, BB=NULL, 
                           columnwise=TRUE, file_genotype=NULL, check=FALSE, workingmemGb=8, csv=FALSE){

 ## check of parameters
 check.inputs(path=path, bin_path=bin_path, file_genotype=file_genotype, workingmemGb=workingmemGb)

 ## checking that AA, AB, and BB have been specified. 
 ## Has AA, AB, BB been assigned numeric values
  if(is.null(AA) |  is.null(BB))
  {
     stop(" AA and BB must be assigned a numeric value \n")
  }

  if(!is.numeric(AA) | !is.numeric(BB))
     stop(" AA and/or BB must be an integer value corresponding to the AA and/or BB genotype, respectively,  in the marker genotype file. ")
  if(!is.null(AB))
     if(!is.numeric(AB))
       stop(" AB must be an integer value corresponding to the AB genotype in the marker genotype file. ")
  ## if there are no hets. 
  if(is.null(AB))
     AB <- 18923282  ## no hets so setting it to something weird



 ## checking for a correct genotype file
 if(check)
   check.genofile(fnameIN=file_genotype, dirPath=path, AA=AA, AB=AB, BB=BB, 
                  check_num_geno_in_row=TRUE, check_genotypes=TRUE, csv=csv)

 if(.Platform$OS.type == "unix") {
    genofile <- paste(path, "/", file_genotype, sep="")
  } else {
   genofile <- paste(path, "\\", file_genotype, sep="")
  }


 if(.Platform$OS.type == "unix") {
    bin_path <- paste(bin_path, "/",  sep="")
  } else {
   bin_path <- paste(bin_path, "\\",  sep="")
  }


  ## Rcpp function to get dimensions of ASCII genotype file
  if (columnwise){
    dim_of_bin_M <- getRowColumn(fname=genofile, csv=csv )
  } else {
    dim_of_bin_M <- getRowColumn(fname=genofile, csv=csv )
    dim_of_bin_M <- c(dim_of_bin_M[2], dim_of_bin_M[1])
  }

  ## Rcpp function to create binrary packed Mt file
#  create.bin.Mt(genofile, bin_path, AA, AB, BB, workingmemGb, dim_of_bin_M, csv)

  ## Rcpp function to create binrary packed M file
#  create.bin.M (genofile, bin_path, AA, AB, BB, workingmemGb, dim_of_bin_M, csv  )


  ## Rcpp function to create binary packed M and Mt file from 
  ## columnwise or non-columnwise data
  create.bin(genofile, bin_path, columnwise, AA, AB, BB, workingmemGb, dim_of_bin_M, csv  )
  
  if(.Platform$OS.type == "unix") {
    binfileM <- paste(path, "/", "M.bin", sep="")
    binfileMt <- paste(path, "/", "Mt.bin", sep="")
  } else {
   binfileM <- paste(path, "\\", "M.bin", sep="")
   binfileMt <- paste(path, "\\", "Mt.bin", sep="")
  }


  return(list("binfileM"=binfileM, "binfileMt"=binfileMt, "dim_of_bin_M" = dim_of_bin_M) )


}




extract_geno <- function(bin_path=NULL, colnum=NULL, workingmemGb=8, dim_of_bin_M=NULL )
  {
    ## Rcpp function to extra a column of genotypes from a  binary packed file of M

   if(.Platform$OS.type == "unix") {
      binfileM <- paste(bin_path, "/", "M.bin", sep="")
    } else {
     binfileM <- paste(bin_path, "\\", "M.bin", sep="")
    }
    selected_locus <- colnum - 1  ## to be consistent with C++'s indexing starting from 0
    geno <- extract_geno_rcpp(f_name_bin=binfileM, max_memory_in_Gbytes = workingmemGb, 
                              selected_locus=selected_locus, dims=dim_of_bin_M)

    return(geno)

  }



constructX <- function(currentX=NULL, loci_indx=NULL, bin_path=NULL, workingmemGb=8, dim_of_bin_M=NULL)
  {
    ## R function to construct the design matrix X
    ## Args
    ##   currentX    current model matrix
    ##   loci        the marker loci to be included as fixed QTL effects (additive model)

   genodat <- extract_geno(bin_path=bin_path, colnum=loci_indx, workingmemGb=workingmemGb, dim_of_bin_M=dim_of_bin_M)

   newX <- cbind(currentX, genodat)
   
   return(newX)

  }


#' @title multiple locus association mapping for genome-wide association studies
#' @description Main  function for performing multiple locus association mapping via whole-genome multi-locus association mapping (WMAM)
#' @param numcores a numeric value for the number of cores that are available for distributed computing. 
#' @param workingmemGb a numeric value. It specifies the amount of memory (in Gigabytes) available for reading analysis. 
#' @param colname.trait  the column name of the trait/response in the phenotypic data file. 
#' &param colname.vars   a character vector containing the column names of 
#'                        the explanatory/independent variables in the phenotypic data file.
#' @param map   the (data frame) object obtained from running \code{\link{read.map}}. If not specifed, a generic map will 
#'              be assumed. 
#' @param pheno  the (data frame) object  obtained  from running \code{\link{read.phenotypes}}.
#' @param geno   the (list) object obtained from running \code{\link{read.genotypes}}.
#' @param alpha  the type 1 error rate where setting \code{alpha} to 0.05 say is a 5\% error rate.
#' @param  error_checking a logical value. When \code{TRUE}, it performs checks of XXXXX of the calculations.  
#'
#' @details
#' The \code{multiple_locus_am} function is an R/Rcpp implementation of whole-genome multi-locus association mapping. The method is a 
#' multiple locus method that is a hybrid of whole-genome and multi-locus association mapping. Briefly, a multiple locus model
#' is built iteatively, by fitting a whole-genome model at each step. It differs from whole-genome mapping methods because we
#' get a parmonious set of marker loci that are in strongest association with a trait.  Also, it differs from multi-locus 
#' association mapping methods because at each iteration of the model building process, we fit all snp loci 
#' simultaneously, as opposed to fitting them one at a time. 
#'
#' NEEDS MORE
#'
#'
#' @seealso \code{\link{read.genotypes}}, and \code{\link{read.phenotypes}}.
#'
#' @return
#' something here .... 
multiple_locus_am <- function(numcores=1,workingmemGb=8, 
                              colname.trait = NULL, 
                              colname.vars  = NULL,
                              map = NULL,
                              pheno=NULL, 
                              geno=NULL, 
                              alpha=0.05, error_checking=FALSE){
 ## Core function for performing whole genome association mapping with EMMA
 ## Args
 ## numcores        number of cores available for computation
 ## memoryGb        maximum amount of working memory available for computation
 ## pheno           data frame or numeric vector. If data frame, then first column is the response and the 
 ##                 remaining columns are explanatory variables to include in the model. If a numeric vector, then it 
 ##                 is only a response to be fitted. 
 ## geno            if geno is a matrix or data frame, then the user has not read.genotypes and a bin packed file
 ##                 has not been created. If it is a character string, then it is the file location of the binary packed files. 
 ## alpha           significance level at which to perform the likelihood ratio test
 ## error_checking  when true, it performs some checks of the calculations


if(.Platform$OS.type == "unix") {
    bin_path  <- paste(bin_path, "/",  sep="")
} else {
   bin_path    <- paste(bin_path, "\\",  sep="")
}


check.inputs(numcores=numcores, workingmemGb=workingmemGb, alpha=alpha)

if(is.null(colname.trait)){
   cat(" Error: the name of the column containing the trait data must be given.")
   stop(" multiple_locus_am has terminated with errors. ")
}

if(is.null(pheno)){
   cat(" Error: the pheno parameter has not been set. \n")
   cat("        Set this parameter to the object that contains \n")
   cat("        the phenotypic data. This object is the result of running  \n")
   cat("        read.phenotypes. \n")
   stop(" multiple_locus_am has terminated with errors.")
}

if(is.null(geno)){
   cat(" Error: the geno parameter has not been set. \n")
   cat("        Set this parameter to the object that contains \n")
   cat("        the phenotypic data. This object is the result of running  \n")
   cat("        read.genotypes. \n")
   stop(" multiple_locus_am has terminated with errors.")
}

## checking list structure of geno
if(!is.list(geno)){
  cat(" Error: the geno object is not a list object. \n")
  cat("       The geno object is obtained from running read.genotypes.\n")
  stop(" multiple_locus_am has terminated with errors.")
}

nms <- names(geno)
indx <- match(nms, c("binfileM", "binfileMt", "dim_of_bin_M"))
if(is.na(indx)){
  cat(" Error: there is a problem with the list structure of the geno object. \n")
  cat("        It should contain the elements binfileM, binfileMt, and dim_of_bin_M. \n")
  stop(" multiple_locus_am has terminated with errors.")
}


if(is.null(map)){
    cat("\n\n  Warning: no map object has been specified. A generic map \n")
    cat("          will be assumed.                                \n\n")
    map <- data.frame(Mrk= paste("M", 1:geno[["dim_of_bin_M"]]), Chrm=1, Pos=1:geno[["dim_of_bin_M"]])
}

   continue <- TRUE
   selected_loci <- NA

   if(is.numeric(pheno)){
     trait <- pheno
     currentX <- matrix(data=1, nrow=length(pheno), ncol=1)  ## column matrix of 1's
   } else {
     ## matrix or data frame
     trait <-  pheno[[1]]
     mf <- paste(names(pheno)[-1], collapse=" + ")
     mf <- paste(" ~ ", mf, sep="")
     mf <- as.formula(mf)
     currentX <- model.matrix(mf, data=pheno)
   }


  while(continue){
       ## Calculate MMt and its inverse
       cat(" --------- Performing iteration of whole-genome search --------------  \n")
       cat(" Calculating MMt \n")
       MMt <- calculateMMt(geno=geno[["binfileM"]], workingmemGb=workingmemGb, numcores=numcores, 
                        dim_of_bin_M = geno[["dim_of_bin_M"]], selected_loci=selected_loci )
       invMMt <- solve(MMt)

      ## perform likelihood ratio test for variance component Var_g
      res_base <- emma.REMLE(y=trait, X= currentX , K=diag(nrow(MMt)) )
      res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100)
      res_full

      ts <- 2 * ( res_full$REML -  res_base$REML )  ## test statistic

      critical_value <- qchisq(1-(2*alpha) , df=1)     ## significance threshold for a mixture of 
                                                       ## chisq distributions of the form 
                                                       ## 0.5 \Chisq_0 + 0.5 \Chisq_1. Here, 
                                                       ## the critical value is found by taking twice
                                                       ## the alpha value under a chi-square distribution with 
                                                       ## 1 degree of freedom. 
   #ts <- 10
     if(ts < critical_value){
        continue <- FALSE
     } else {  ## QTL present
        cat(" Calculating H \n")
        H <- calculateH(MMt=MMt, varE=res_full$ve, varG=res_full$vg)
        cat(" Calculating P \n")
        P <- calculateP(H=H, X=currentX)
        cat(" Calculating MMt_sqrt_and_sqrtinv \n")
        MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking, verbose=FALSE)

    #    hat_a <- calculate_reduced_a(varG=res_full$vg, bin_path=bin_path, P=P, y=trait, 
    #                                 workingmemGb=workingmemGb, dim_of_bin_M=geno[["dim_of_bin_M"]], 
    #                                 selected_loci=selected_loci)
         ## hat_a <- calculate_reduced_a(varG=res_full$vg, P=P, MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], y=trait)
         cat(" calculating h_a \n")
         hat_a <- calculate_reduced_a(varG=res_full$vg, P=P, MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], y=trait)  

        cat(" var_hat_a \n")
        var_hat_a    <- calculate_reduced_vara(X=currentX, varE=res_full$ve, varG=res_full$vg, invMMt=invMMt, 
                                            MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]])
   
        cat(" a_and_vara \n")
        a_and_vara  <- calculate_a_and_vara(bin_path=bin_path,  maxmemGb=workingmemGb, dims=geno[["dim_of_bin_M"]],
                         selectedloci = selected_loci,
                         invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]], transformed_a=hat_a, transformed_vara=var_hat_a)
  
        ## outlier test statistic
        tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
        indx <- which(tsq == max(tsq))   ## index of largest test statistic. However, need to account for other loci 
                                         ## already having been removed from M which affects the indexing

        ## taking first found qtl
        indx <- indx[1]

        orig_indx <- seq(1, geno[["dim_of_bin_M"]][2])  ## 1:ncols
        cat(" indx = ", indx, "\n")
        new_selected_locus <- orig_indx[indx]
        print(new_selected_locus)
        selected_loci <- c(selected_loci, new_selected_locus)
        cat(" Selected loci = ", selected_loci, " \n")
        currentX <- constructX(currentX=currentX, loci_indx=new_selected_locus, bin_path = bin_path, 
                                      dim_of_bin_M=geno[["dim_of_bin_M"]]) 
        cat("   ....  Found qtl at index ... ", new_selected_locus, " \n")
     }  ## if ts < critical_value



  }  ## end while continue




return(selected_loci)

} ## end multiple_locus_am















