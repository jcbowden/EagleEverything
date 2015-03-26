### To run multiple CPU
## salloc --ntasks-per-node=1  --ntasks-per-node=6 --mem=20gb --time=30:0 srun --pty bash

## To DO
## how to handle header in genotype file
## how to handle missing genotypes - are these allowed? probably not. 
## tidy up emma code
## how to write the model matrix when pheno can contain an arbritary number of columns

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
#' @param check_num_geno_in_row a logical value. When \code{TRUE}, the number of genotypes in each row is returned.
#' @param check_genotypes a logical value. When \code{TRUE}, it checks that the marker genoyptes are either 0, 1, or 2. 
#'
#' @return  a list is returned with elements \code{file_exists} and \code{\num_genotypes_per_row}. \code{num_genotypes_per_row}
#'          will be \code{NULL} unless the \check_num_geno_in_row} parameter has been set to \code{TRUE}.
#' @seealso \code{\link{read.genotypes}}
check.genofile <- function(fnameIN=NULL, dirPath=getwd(), 
                           check_num_geno_in_row=FALSE,
                           check_genotypes=FALSE)
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


  ## has deletion file name been provided
  if(is.null(fnameIN))
      stop(" No name of the genotype file has been supplied.")


  ## does the file contain a constant number of records. 
  if(check_num_geno_in_row){
     res[["num_genotypes_per_row"]] <- count.fields(file=file_geno)
     if(length(unique( res[["num_genotypes_per_row"]] ))>1){
       cat(" File ", file_geno, " has  differing numbers of genotypes per line ... \n")
       cat(" The differing number of gentoypes per line are ", unique(cf), "\n")
       stop(" Please modify the genotype file to have the same number of genotypes per row. \n")
     } else {
       cat(" File ", file_geno, " has ", unique( res[["num_genotypes_per_row"]] ), " genotypes per line...\n\n")
     } ## end if else
  }  ## end if checkrownum


  ## does the file contain 0,1,2 genotypes only (implemented in Rcpp for speed)
  if(check_genotypes)
      checkGenotypes(file_geno)

  return(res)
}






calculateMMt <- function(geno=NULL, workingmemGb, numcores, selected_loci=NA, dim_of_M=NULL)
{
 ## R interface to Rcpp code to calculate M %*% t(M)
 ## Args
 ##      geno        absolute path + file name of binary packed M file
 ##      workingmemGb    amount of memory in Gbytes available for creation of MMt
 ##      numcores    number of cores for matrix operations
 ##      selectedloci an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                   have been selected to act as fixed QTL effects in the model. 
 ##      dim_of_M    numeric vector with the row, column numbers of M. 
  #------------------------------------------
  # bin file about to be overwritten
  #------------------------------------------

  if(!file.exists(geno)){
    cat(" Error: The binary packed file ", geno, " cannot be found.\n")
    stop(" calculateMMt has terminated with errors.") 
   }


  MMt <- calculateMMt_rcpp( f_name_bin=geno, selected_loci = selected_loci,
                               max_memory_in_Gbytes=workingmemGb, num_cores=numcores, 
                               dims= dim_of_M)
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



mistake_calculate_reduced_a <- function(varG=NULL, bin_path=getwd(), P=NULL, y=NULL, workingmemGb=8, dim_of_M=NULL, 
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
                                 dims=dim_of_M , selected_loci = selected_loci )





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
#' @details
#' A space separated ASCIII file is assumed. This file is allowed to contain \code{NA} values (coded as \code{NA}). However, if
#' \code{NA} values are contained in the covariates (or explanatory variables) and these covariates are being used by 
#' \code{\link{multiple_locus_am}}, then an error will occur. \code{NA} values are allowed in the trait values.  
#' @seealso \code{\link{read.genotypes}}
#' @return 
#' a data frame is returned of the phenotypic data. If \code{header} is true, the 
#' names of the columns will be as specified by the first row of the phenotypic file. If \code{header} is \code{FALSE}, 
#' generic names are supplied by R in the form of V1, V2, etc.   
#'
read.phenotypes <- function(path=getwd() , file_phenotype = NULL, header=TRUE){

  check.inputs(path=path, file_phenotype=file_phenotype)

  if(.Platform$OS.type == "unix") {
    phenofile <- paste(path, "/", file_phenotype, sep="")
  } else {
   phenofile <- paste(path, "\\", file_phenotype, sep="")
  }

  phenos <- read.table(phenofile, header=header)

  ## produce summary of phenotypic information 
  cat("\n\n", sprintf("     Summary of contents of phenotypic file "), "\n")
  cat(sprintf("  ---------------------------------------------------------- "), "\n")
  cat(sprintf("  File name of phenotypic file is : %s", phenofile), "\n")
  cat(sprintf("  File has been read in as a data frame"), "\n")
  cat(sprintf("  Number of columns read: %20d", ncol(phenos)), "\n")
  cat(sprintf("  Number of rows read:    %20d", nrow(phenos)), "\n")
  if(header){
    cnames <- paste(names(phenos), collapse="   ")
    cat(sprintf("  Column names are :    %s", cnames), "\n")
  }
  cat(sprintf("  Column 1 is the response and is called %s", names(phenos)[1]), "\n")
  if(ncol(phenos) > 1)
  {
     for(ii in 2:ncol(phenos))
     {
        cat(sprintf("  Column %d is an explanatory variable and is called   %s   and is of class %s ", ii, names(phenos)[ii], class(phenos[[ii]])), "\n")
     } 
  }
  cat("\n", sprintf( " Warning!!!  If the classes of these explanatory variables is not correct, these will need to be changed by the user."), "\n\n")

  return(phenos)


}





create.bin.Mt <- function(file_genotype, bin_path, workingmemGb, dim_of_M){
 ## an Rcpp function to create the packed binary file of the genotype data Mt
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## bin_path         path name for where binary files are to be stored. 
 ## workingmemGb     available memory for converstion to packed binary
 ## dim_of_M             row, column dimensions of M.  

 binMtfile <- paste(bin_path, "Mt.bin", sep="") ## file name for binary packed Mt file.

 createMt_rcpp(f_name = file_genotype, f_name_bin = binMtfile,  
               max_memory_in_Gbytes=workingmemGb,  dims = dim_of_M )

 cat(" A packed binary file called ", binMtfile, " has been created for the transform of the genotype data.\n")


}



create.bin.M <- function(file_genotype, bin_path, workingmemGb, dim_of_M){
 ## an Rcpp function to create the packed binary file of the genotype data M
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## bin_path         path name for where binary files are to be stored. 
 ## workingmemGb     available memory for converstion to packed binary
 ## dim_of_M             row, column dimensions of M.  

 binMfile <- paste(bin_path, "M.bin", sep="") ## file name for binary packed Mt file.

 createM_rcpp(f_name = file_genotype, f_name_bin = binMfile,
               max_memory_in_Gbytes=workingmemGb,  dims = dim_of_M )

 cat(" A packed binary file called ", binMfile, " has been created for the genotype data.\n")



}



#' @title Read  marker genotype file.
#' 
#' @description
#' \code{read.genotypes} reads in the genotypes from the marker genoytpe file.
#' @param path  character vector containing the directory path to where the marker genotype file is located.
#' @param bin_path  character vector containing the directory path to where the binary converted marker genotype files are to be located. 
#' @param file_genotype character vector containing the name of the marker genotype file.
#' @param  check  a logical value. If \code{TRUE}, then \code{\link{check.genofile}} is run to check the marker genotype file
#'         for errors. 
#' @param workingmemGb a numeric value. It specifies the amount of memory (in Gigabytes) available for reading in the marker 
#'                     genotype file. 
#' @details
#' \code{read.genotypes} reads the marker genoytpes from file and converts it into two binary files:
#' a binary packed file of the original marker data and a binary packed file of the transpose of the original 
#' marker data.  
#'
#' The two binary files created by \code{read.genotypes} are called M.bin and Mt.bin and they are created in 
#' directory \code{bin_path}. 
#'
#' It is assumed that the marker genotype file is a space separated ASCII file with the rows being the samples/individuals 
#' and the columns being the SNP loci. 
#'
#' Missing values are not allowed. If present, this will cause \code{read.genotypes} to fail. Such values should be imputed to 
#'      0,1, or 2 values first. 

#' The \code{path} and \code{bin_path} parameters, if not set, will default to the current directory. 
#'
#' The \code{workingmemGb} parameter should be set to the largest amout of RAM for the machine, when practical.  If the marker 
#' genotype file is too large for the available memory, then it converted, chunk-wise,  into a single binary packed file.  
#' The consequences of reading the marker genotype file chunk- or block-wise is a longer processing time. 
#'
#' @return  a list is returned with elements \code{binfileM} , \code{binfileMt}, and \code{dim_of_M}. 
#' which is the absolute path and file name of the binrary file for the marker genotype data, 
#' the absolute path and file name of the binary file for the transpose of the marker genotype data, 
#' and the number of rows and columns in the marker genotype file, respectively.
#' @seealso \code{\link{read.genofile}}
read.genotypes <- function(path=getwd(), bin_path=getwd(), file_genotype=NULL, check=FALSE, workingmemGb=8){
 ## an Rcpp based function to check for problems in reading the genotyope file 
 ## and to created binary packed M and Mt. 

 ## bin_path      directory in which binary files are to be created. Default is the working directory
 ## path          directory containing raw ASCII genotype file. 
 ## file_genotype file name of the genotype file
 ## check         if true, then some simple checks of the genotype file are performed such as the number of 
 ##               genotypes per row and if values different from 0,1,2 are present. 
 ## header        FUNCTIONALITY NOT WORKING YET - NOT SURE HOW TO DEAL WITH THIS
 ## workingmemGb  amount of memory available for converting data to binary packed form.
 #3
 ## Returns
 ##  The absolute path and file name of the binary packed genotype file
 ## check of parameters
 check.inputs(path=path, bin_path=bin_path, file_genotype=file_genotype, workingmemGb=workingmemGb)

 ## checking for a correct genotype file
 if(check)
   check.genofile(fnameIN=file_genotype, dirPath=path, check_num_geno_in_row=TRUE, check_genotypes=TRUE)

 if(.Platform$OS.type == "unix") {
    genofile <- paste(path, "/", file_genotype, sep="")
  } else {
   genofile <- paste(path, "\\", file_genotype, sep="")
  }

  ## Rcpp function to get dimensions of ASCII genotype file
  dim_of_M <- getRowColumn(fname=genofile)


  ## Rcpp function to create binrary packed Mt file
  create.bin.Mt(genofile, bin_path, workingmemGb, dim_of_M)

  ## Rcpp function to create binrary packed M file
  create.bin.M (genofile, bin_path, workingmemGb, dim_of_M  )

  
  if(.Platform$OS.type == "unix") {
    binfileM <- paste(path, "/", "M.bin", sep="")
    binfileMt <- paste(path, "/", "Mt.bin", sep="")
  } else {
   binfileM <- paste(path, "\\", "M.bin", sep="")
   binfileMt <- paste(path, "\\", "Mt.bin", sep="")
  }


  return(list("binfileM"=binfileM, "binfileMt"=binfileMt, "dim_of_M" = dim_of_M) )


}




extract_geno <- function(bin_path=NULL, colnum=NULL, workingmemGb=8, dim_of_M=NULL )
  {
    ## Rcpp function to extra a column of genotypes from a  binary packed file of M

   if(.Platform$OS.type == "unix") {
      binfileM <- paste(bin_path, "/", "M.bin", sep="")
    } else {
     binfileM <- paste(bin_path, "\\", "M.bin", sep="")
    }
    selected_locus <- colnum - 1  ## to be consistent with C++'s indexing starting from 0
    geno <- extract_geno_rcpp(f_name_bin=binfileM, max_memory_in_Gbytes = workingmemGb, 
                              selected_locus=selected_locus, dims=dim_of_M)

    return(geno)

  }



constructX <- function(currentX=NULL, loci_indx=NULL, bin_path=NULL, workingmemGb=8, dim_of_M=NULL)
  {
    ## R function to construct the design matrix X
    ## Args
    ##   currentX    current model matrix
    ##   loci        the marker loci to be included as fixed QTL effects (additive model)

   genodat <- extract_geno(bin_path=bin_path, colnum=loci_indx, workingmemGb=workingmemGb, dim_of_M=dim_of_M)

   newX <- cbind(currentX, genodat)
   
   return(newX)

  }


#' @title multiple locus association mapping for genome-wide association studies
#' @description Main  function for performing multiple locus association mapping via whole-genome multi-locus association mapping (WMAM)
#' @param numcores a numeric value for the number of cores that are available for distributed computing. 
#' @param workingmemGb a numeric value. It specifies the amount of memory (in Gigabytes) available for reading analysis. 
#' @param bin_path  character vector containing the directory path to where the binary converted marker genotype files are located. 
#' @param geno the name of the list object returned from running \code{\link{read.genotypes}}.
#' @param pheno  the name of the data frame object  returned  from running \code{\link{read.phenotypes}}.
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
multiple_locus_am <- function(numcores=1,workingmemGb=8, bin_path=getwd(), geno=NULL,   pheno=NULL, alpha=0.05, error_checking=FALSE){
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


check.inputs(numcores=numcores, workingmemGb=workingmemGb, bin_path=bin_path, alpha=alpha)

if(is.null(geno)){
  cat(" Error: The geno parameter must be set to the output from read.genotypes(). \n")
  stop(" multiple_locus_am has terminated with errors. ")
} else {
  if(!is.list(geno)){
   cat(" Error: geno parameter must be a list object created from the output of read.genotypes(). \n")
   stop(" multiple_locus_am has terminated with errors.")
  }  else {
    if(!file.exists(geno[["binfileM"]])){
      cat(" Error: the binary packed file ", geno[["binfileM"]]," could not be found. \n")
      cat("        This file is created by running read.genotypes(), saving its \n")
      cat("        contents to an object, and setting the geno parameter to this object. \n")
      stop(" multiple_locus_am has terminated with errors.")
     }
  }
}

if(is.null(pheno)){
   cat(" Error: the pheno parameter has not been set. \n")
   cat("        Set this parameter to the object that contains \n")
   cat("        the phenotypic data. This object should be of class \n")
   cat("        data.frame, matrix, or numeric. \n")
   stop(" multiple_locus_am has terminated with errors.")
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
                        dim_of_M = geno[["dim_of_M"]], selected_loci=selected_loci )
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
    #                                 workingmemGb=workingmemGb, dim_of_M=geno[["dim_of_M"]], 
    #                                 selected_loci=selected_loci)
         ## hat_a <- calculate_reduced_a(varG=res_full$vg, P=P, MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], y=trait)
         cat(" calculating h_a \n")
         hat_a <- calculate_reduced_a(varG=res_full$vg, P=P, MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], y=trait)  

        cat(" var_hat_a \n")
        var_hat_a    <- calculate_reduced_vara(X=currentX, varE=res_full$ve, varG=res_full$vg, invMMt=invMMt, 
                                            MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]])
   
        cat(" a_and_vara \n")
        a_and_vara  <- calculate_a_and_vara(bin_path=bin_path,  maxmemGb=workingmemGb, dims=geno[["dim_of_M"]],
                         selectedloci = selected_loci,
                         invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]], transformed_a=hat_a, transformed_vara=var_hat_a)
  
        ## outlier test statistic
        tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
        indx <- which(tsq == max(tsq))   ## index of largest test statistic. However, need to account for other loci 
                                         ## already having been removed from M which affects the indexing

        ## taking first found qtl
        indx <- indx[1]

        orig_indx <- seq(1, geno[["dim_of_M"]][2])  ## 1:ncols
        cat(" indx = ", indx, "\n")
        new_selected_locus <- orig_indx[indx]
        print(new_selected_locus)
        selected_loci <- c(selected_loci, new_selected_locus)
        cat(" Selected loci = ", selected_loci, " \n")
        currentX <- constructX(currentX=currentX, loci_indx=new_selected_locus, bin_path = bin_path, 
                                      dim_of_M=geno[["dim_of_M"]]) 
        cat("   ....  Found qtl at index ... ", new_selected_locus, " \n")
     }  ## if ts < critical_value



  }  ## end while continue




return(selected_loci)

} ## end multiple_locus_am















