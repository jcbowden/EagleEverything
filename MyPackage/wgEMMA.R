# Warranty
# This software is distributed under the GNU General Public License.
#
#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##              Checks 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##  Phenotype data
##     differing number of rows in pheno file and geno file.          x
##     phenotype file not found                                       x
##     checking if code runs with NAs in pheno                        x


## Genotype data
##     differing number of loci in map file and geno file.            x
##     differing number of elements in a row for the genotype file    x
##     genotype file not found                                        x
##     genotype text file has funny characters                        x
##     csv = TRUE                                                     x

## Map file
##     map file not found                                             x

## AM
##      not set pheno                                                 x
##      not set geno                                                  x
##      not set map                                                   x
##      trait name incorrect                                          x
##      trait name not given                                          x
##      fixed effect names wrong                                      x
##      pheno object misspecified                                     x
##      test if objects supplied to AM exist                          x
##      geno object not a list object                                 x
##      pheno object not a data frame




#' Package documentation
#'
#' @name AMplus-package
#' @title A package for genome-wide association mapping with multiple-locus models and for populations of 
#' arbitrary structure. 
#' @docType package
#' @author Andrew W. George \email{andrew.george@@csiro.au}
#' @details See below
#' @description A fully documented R package for performing association mapping with 
#' multiple-locus models on a genome-wide scale. 
#' @section Motivation:  Data from genome-wide association studies are analysed, commonly, with single 
#' locus models, that is,  on a locus-by-locus basis. Association mapping with multiple-locus 
#' models is more powerful but most multiple-locus methods do not scale well with study size
#' and can be difficult to implement in practice.  This package was designed to make 
#' genome-wide association mapping with multiple-locus models simple and practical. 
#' It is based on linear mixed models so is best suited for 
#' traits which are continuous/ quantitative. However, 
#' association mapping with binary traits is also possible.
#'
#'@section Reading in  Marker Data: Our package can accept genotype data in 
#' plain text files and allelic data in PLINK ped files. These files can be 
#' larger than the memory capacity of the machine. Other formats 
#' can also be handled by using plink (with the \code{recode} option) to 
#' convert them into ped files.  See \code{\link{ReadMarker}} for more details.
#'
#' @section Quick start guide:
#' \itemize{
#' \item{}{Create a marker file (plain text or PLINK ped format) with no column 
#' headings and where the rows are the individuals and the columns are the marker loci.}
#' \item{}{Create a map file (optional) that has three columns, the name of the 
#' marker locus, the chromosome number, and the map location of the marker locus. 
#' Column headings are required.}
#' \item{}{Create a phenotype file with the trait data 
#'    and if appropriate, explanatory variable data. Column headings are required. }
#' \item{}{Input marker data with \code{\link{ReadMarker}} }
#' \item{}{Input the map file, if present, with \code{\link{ReadMap}}}
#' \item{}{Input the phenotypic data with \code{\link{ReadPheno}}}
#' \item{}{Perform genome wide association mapping with \code{\link{AM}}}    
#' }
#'
#' @section Output: The aim of a GWAS is to identify those marker loci 
#' closest to the genes that are influencing the trait. So, when the GWAS data are 
#' analysed, a set of marker loci labels are returned as the output. These marker 
#'  loci are closest to the genes underlying the trait and are found while simultaneously 
#' accounting for multiple marker-trait associations, familial relatedness, and 
#' fixed effects such as population structure (if included in the analysis).
#' More detailed output such as the additive effect of the marker locus, its 
#' significance in the multiple-locus model ( measured by a p-value), and 
#' an estimate of the amount of variation explained by the locus can be 
#' obtained by running the summary function (\code{\link{SummaryAM})}.
#'
#'@section Where to get help: To see an overview of the package and its functions type
#' \preformatted{
#' library(, AMplus)
#' }
#'
#' For detailed help on a function, type
#'
#' \preformatted{
#' help(functionname)
#' }
#'
#' where \emph{functionname}  is the name of the function.
#'
#' A YouTube quick start presentation can be found here ... 
#'
#' Powerpoint slides found here Website
#'
#' @keywords Association mapping, multiple-locus models, linear mixed models.
NULL









## internal function to get absolute path name under unix and windows
## but only if absolute path has not already been specified. 
fullpath <- function(fname){
 ## check if full path has been given
  filen <- fname  ## initialize
  if (! (length(grep("/", fname)) > 0 || length(grep("[\\]", fname)) > 0 ) ){
     if(.Platform$OS.type == "unix") {
       filen <- paste(getwd(), "/", fname, sep="")
     } else {
      filen <- paste(getwd(), "\\", fname, sep="")
     }
  }
  return(filen)
}


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

 emma.eigen.L.w.Z <- function (Z, K, complete = TRUE, ngpu=0) 
{
    if (complete == FALSE) {
        vids <- colSums(Z) > 0
        Z <- Z[, vids]
        K <- K[vids, vids]
    }
    res <- K %*% crossprod(Z, Z)
    ## cannot use eigen_mgpu here because matrix is not symmetric
    eig <- eigen(res, symmetric = FALSE, EISPACK = TRUE)
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
    esp = 1e-10, eig.L = NULL, eig.R = NULL, ngpu=0) 
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
            eig.L <- emma.eigen.L.wo.Z(K, ngpu)
        }
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.wo.Z(K, X, ngpu)
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
            eig.L <- emma.eigen.L.w.Z(Z, K, ngpu)
        }
        if (is.null(eig.R)) {
            eig.R <- emma.eigen.R.w.Z(Z, K, X, ngpu)
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





emma.eigen.L.wo.Z <- function (K, ngpu=0) 
{  
    if(ngpu > 0){
      if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
         eig <- rcppMagmaSYEVD::eigen_mgpu(K, symmetric=TRUE)
       }

     } else {
      eig <- eigen(K, symmetric = TRUE)
     }
    return(list(values = eig$values, vectors = eig$vectors))
}

emma.eigen.R.wo.Z <-  function (K, X, ngpu=0) 
{
    n <- nrow(X)
    q <- ncol(X)
    dn <- diag(n)
    S <- dn - X %*% solve(crossprod(X, X)) %*% t(X)
    gc()
    if(ngpu > 0){
     if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
       eig <- rcppMagmaSYEVD::eigen_mgpu(S %*% (K + dn) %*% S, symmetric = TRUE, only_values=FALSE)
     }
    } else {
       eig <- eigen(S %*% (K + dn) %*% S, symmetric = TRUE)
    }


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



 emma.REMLE <-  function (y, X, K, Z = NULL, ngrids = 100, llim = -10, ulim = 10,  ngpu=0,
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
            eig.R <- emma.eigen.R.wo.Z(K, X, ngpu)
        }
        etas <- crossprod(eig.R$vectors, y)
        logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
        m <- length(logdelta)
        delta <- exp(logdelta)
#print("eig.R$values")
#print(eig.R$values)
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
 #print(c(i, dLL[i], dLL[i+1], esp, m-1))
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


check.for.NA.in.trait <- function(trait=NULL)
{
     ## internal function for AM 
     ## to return the positions of NA in a trait

       ## check for NA's in trait
        indxNA <- which(is.na(trait))
        if(length(indxNA)==0){
          indxNA <- vector("numeric", 0)
        } else {
          ## place in reverse order
          indxNA <- sort(indxNA, decreasing = TRUE)
cat("\n\n WARNING!!!! The individuals in rows ", indxNA, " either have missing trait data \n")
cat("             and/or missing explanatory variable values. These individuals have \n")
cat("             been removed from the analysis.  \n\n")
          if(any(is.na(indxNA))){
            cat("Error:  (internal).  indxNA contains NA values. \n")
            stop(" AM has terminated with errors. ", call. = FALSE)
          }
        }

      return(indxNA)
   } 


check.inputs.mlam <- function (ncpu, availmemGb, colname.trait, colname.feffects, map, pheno, 
                  geno )
{



if(is.null(colname.trait)){
   cat("Error: the name of the column containing the trait data must be given. \n")
   return(TRUE)
}

if(is.null(pheno)){
   cat("Error: the pheno parameter has not been set. \n")
   cat("       Set this parameter to the object that contains \n")
   cat("       the phenotypic data. This object is the result of running  \n")
   cat("       ReadPheno. \n")
   return(TRUE)
}

if(is.null(geno)){
   cat("Error: the geno parameter has not been set. \n")
   cat("       Set this parameter to the object that contains \n")
   cat("       the phenotypic data. This object is the result of running  \n")
   cat("       ReadMarker. \n")
   return(TRUE)
}


if(class(try(class(geno), silent=TRUE)) == "try-error"){
   cat("Error: the object supplied to the geno parameter does not exist. \n")
   cat("       This object is set by running ReadMarker. Type help(ReadMarker) for help \n")
   cat("       on running this command. \n")
   return(TRUE)
}


if(class(try(class(pheno), silent=TRUE)) == "try-error"){
   cat("Error: the object supplied to the pheno parameter does not exist. \n")
   cat("       This object is set by running ReadPheno. Type help(ReadPheno) for help. \n")
   cat("       on running this command. \n")
   return(TRUE)
}


## checking list structure of geno
if(!is.list(geno)){
  cat("Error: the geno object is not a list object. \n")
  cat("     The geno object is obtained from running ReadMarker.Type help(ReadMarker) for help. \n")
  return(TRUE)
}

## checking if pheno is a data frame 
if(!is.data.frame(pheno)){
  cat("Error: the pheno object is not a data frame. \n")
  cat("      It is a ", class(pheno), "\n")
  return(TRUE)
}



nms <- names(geno)
indx <- match(nms, c("binfileM", "binfileMt", "dim_of_bin_M" ))
if(any(is.na(indx))){
  cat("Error: there is a problem with the list structure of the geno object. \n")
  cat("       It should contain the elements binfileM, binfileMt, and dim_of_bin_M. \n")
  cat("       The object supplied contains the elements ", names(geno), "\n")
  return(TRUE)
}

if(is.null(map)){
    cat("WARNING: no map object has been specified. A generic map \n")
    cat("         will be assumed.                                \n\n")
    map <- data.frame(Mrk= paste("M", 1:geno[["dim_of_bin_M"]][2]), Chrm=1, Pos=1:geno[["dim_of_bin_M"]][2])
}






 ## checks for colname.trait
 if(is.null(colname.trait)){
    cat("Error: the column name for the trait/response has not been specified.\n")
    cat("       Please set trait to the column name of the trait data in \n")
    cat(" .     the phenotypic file. The allowable column names are ", names(pheno), "\n")
    return(TRUE)
 }

 if(length(colname.trait)>1){
    cat("Error: multiple column names for the trait have been specified. \n")
    cat("       Only a single column name should be  assigned to trait. \n")
    return(TRUE)
 }

 indx <- match(colname.trait, names(pheno))
 if(any(is.na(indx))){
   cat("Error: the trait column name does not match any of the column names in the phenotypic file. \n")
   cat("       The name that has been supplied is ", pheno, "\n")
   cat("       The column names of the phenotypic file are ", names(pheno), "\n")
   return(TRUE)
 }


 ## checks for colname.feffects
 if(is.null(colname.feffects)){
    cat("\n\n WARNING: no fixed effects have been specified. \n\n")
 }


 indx <- match(colname.feffects, names(pheno))
 if(any(is.na(indx))){
   cat("Error: the feffects option contains column names that do not \n")
   cat("       match any of the column names in the phenotypic file. \n")
   cat("       The column names in the phenotypic file are: \n")
   cat(c("       ", names(pheno), "\n"))
   return(TRUE)
 }



 ## check that geno and pheno contain the same number of individuals
 if(nrow(pheno) !=  geno[["dim_of_bin_M"]][1])
 {
   cat("Error: the number of individuals specified in the phenotypic file is ", nrow(pheno),  "\n")
   cat("       the number of individuals specified in the genotypic file is ",  geno[["dim_of_bin_M"]][1], "\n")
   cat("       The number of individuals should be the same in the two files.\n")
   return(TRUE)
 }

 ## check that map and geno contain the same number of snp
 if(nrow(map) != geno[["dim_of_bin_M"]][2])
 {
   cat("Error: the number of marker loci in the map file is ", nrow(map), "\n")
   cat("       The number of marker loci in the genotypic file is ", geno[["dim_of_bin_M"]][2], "\n")
   cat("       The number of marker loci in the two files should be the same. \n")
   return(TRUE)
 }



  return(FALSE)

}







#### To run multple GPU's
### > export OMP_NUM_THREADS=$PBS_NUM_PPN
## > module load cuda/6.0 R
## > module load R/3.0.0
## > LD_PRELOAD=libnvblas.so R  
## monitoring gpu usage
##    nvidia-smi -l 3


##  library('Rcpp')
##  library('RcppEigen')
##  library('matrixcalc')
##  library('Matrix')
## This builds a dll for the function
## sourceCpp("/home/geo047/gitHUB_WMAM/MyPackage/RcppFunctions.cpp", rebuild=TRUE, quiet=TRUE)
## source("/home/geo047/gitHUB_WMAM/MyPackage/wgEMMA.R")
## source("/home/geo047/gitHUB_WMAM/MyPackage/multiple_am.R")


##-------------------------------
## R Function Declaration 
##-------------------------------


calculateMMt <- function(geno=NULL, availmemGb, ncpu, selected_loci=NA, dim_of_bin_M=NULL, quiet = 0)
{
 ## R interface to Rcpp code to calculate M %*% t(M)
 ## Args
 ##      geno        absolute path + file name of binary packed M file
 ##      availmemGb    amount of memory in Gbytes available for creation of MMt
 ##      ncpu    number of cores for matrix operations
 ##      selectedloci an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                   have been selected to act as fixed QTL effects in the model. 
 ##      dim_of_bin_M    numeric vector with the row, column numbers of M. 
  #------------------------------------------
  # bin file about to be overwritten
  #------------------------------------------


  if(!file.exists(geno)){
    cat(" Error: The binary packed file ", geno, " cannot be found.\n")
    stop(" calculateMMt has terminated with errors.", call. = FALSE) 
   }
  if(!any(is.na(selected_loci))) selected_loci <- selected_loci-1
  MMt <- calculateMMt_rcpp( f_name_bin=geno, selected_loci = selected_loci,
                               max_memory_in_Gbytes=availmemGb, num_cores=ncpu, 
                               dims= dim_of_bin_M, quiet = quiet) 
  return(MMt)

}  ## end function









calculateMMt_sqrt_and_sqrtinv <- function(MMt=NULL, checkres=TRUE, 
                                           quiet = 0 , ngpu=0)
{
  ## R function for calculating the square root of M * M^t
  ## and the inverse of the square root of MMt
  ## where M * M^t has already been created. 
  ## Using SVD for calculation
  ## Args
  ##  MMt   a matrix of M * M^t
  ##  checkres  when true, the accuracy of the inversion is checked. 

  ## testing that MMt is postive definite
  if(!is.positive.definite(MMt)){
    cat(" Error: the matrix multiplication M %*% t(M) is not positive definite. \n")
    cat("        This can occur if there are individuals with identical marker \n")
    cat("        information. Please remove individuals with indentical marker \n")
    cat("        information, remembering also to remove their associated phenotypic \n")
    cat("        information as well. \n")
    stop(" Internal function: calculateMMt_sqrt_and_sqrtinv has terminated with errors.\n", call. = FALSE)
  } 
   res <- list()

   if(ngpu == 0){
      MMt.eigen <- eigen(MMt, symmetric=TRUE )
      sqrt_evals <- diag(sqrt(MMt.eigen$values))
      res[["sqrt"]] <- MMt.eigen$vectors %*% sqrt_evals %*% t(MMt.eigen$vectors)
      rm(MMt.eigen, sqrt_evals)
      gc()
      res[["invsqrt"]] <- chol2inv(chol(res[["sqrt"]]))
   }  else {
      #res <- rcppMagmaSYEVD::sqrt_invsqrt(MMt, symmetric=TRUE)
      if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
        res <- rcppMagmaSYEVD::sqrt_invsqrt(MMt, symmetric=TRUE)
      }
   } 



   if(checkres){
       a <- (res[["sqrt"]] %*% res[["invsqrt"]] )
       if(trunc(sum(diag(a))) != nrow(MMt))
       {
         cat(" \n\n\nWARNING: these results may be unstable.\n")
         cat(" The sum of the diagonal elements of the square root of M %*% t(M) and its inverse is ", sum(diag(a)), " where \n")
         cat("  it should have been ", nrow(MMt), "\n")
         cat("  This can occur if the genotype file contains near identical rows and/or columns.  Please check.\n\n")
      

       } 
   }   ## end if(checkres)
   res <- list(sqrt_MMt=res[["sqrt"]], inverse_sqrt_MMt=res[["invsqrt"]] )
  


} ## end function



calculateH <- function(MMt=NULL, varE=NULL, varG=NULL )
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
    stop(" The varE (residual variance) must be numeric.", call. = FALSE)

  if(varE < 0)
    stop(" VarE cannot be negative.", call. = FALSE)

  if(!is.numeric(varG))
    stop(" The varG (genotypic variance) must be numeric.", call. = FALSE)

  if(varG < 0)
    stop(" VarG cannot be negative.", call. = FALSE)


  if(is.null(MMt))
    stop("MMt cannot be null.")

  return( varE * diag(nrow(MMt)) + varG * MMt)


}


calculateP  <- function(H=NULL, X=NULL, ngpu=0)
{
  ## R function to calculate P matrix
  ## Args:
  ##       H is the variance matrix
  ##       X is the design matrix supplied by the user
  ##       ngpu for the number of gpu
  ## Returns:
  ##   matrix object P

  if(is.null(H))
    stop(" H must be specified.", call. = FALSE)
  if(is.null(X))
    stop(" A design matrix has not be specified. ", call. = FALSE)

   if(nrow(H) != nrow(X))
      stop(" The number of rows in H and X are not the same.", call. = FALSE)

 Hinv <- chol2inv(chol(H))
 P <- Hinv - Hinv %*% X %*% solve( t(X) %*% Hinv %*% X )  %*% t(X) %*% Hinv

  return(P)

}


calculate_reduced_a <- function(varG=NULL, P=NULL, MMtsqrt=NULL, y=NULL, quiet=FALSE)
{

  if( !(nrow(P) ==  length(y))){
    cat(" Error:  there is a problem with the  dimensions of  P, and/or the vector y.")
    cat("         They should  be of the dimension (n x n), and a vector of length n.")
    cat(" The dimensions are: \n")
    cat(" dim(P)      = ", dim(P), "\n")
    cat(" length(y)   = ", length(y), "\n")
    stop(call. = FALSE)

  }

 if(is.null(varG))
   stop(" VarG must be specified.", call. = FALSE)

  if(is.null(P))
   stop(" P must be specified", call. = FALSE)


  if(is.null(y))
   stop(" y must be specified", call. = FALSE)

    a <- varG * MMtsqrt %*% P %*% y

return(a)

}



mistake_calculate_reduced_a <- function(varG=NULL, P=NULL, y=NULL, availmemGb=8, dim_of_bin_M=NULL, 
                                 selected_loci=NA, quiet = 0)
{
 ## Rcpp function to calculate the BLUP (a) values under a dimension reduced model
 ## Args:
 ##    varG is a scalar value
 ##    P   is a n x n matrix
 ##    y   is a n x 1 vector
 ##
 ## a* = sigma^2_a * t(M) * P * y
 ## Returns:
 ##   a numeric vector of dimension reduced a values 

 if(is.null(varG))
   stop(" VarG must be specified.", call. = FALSE)

  if(is.null(P))
   stop(" P must be specified", call. = FALSE)

 
  if(is.null(y))
   stop(" y must be specified", call. = FALSE)

 
  if( !(nrow(P) ==  length(y))){
    cat(" Error:  there is a problem with the  dimensions of  P, and/or the vector y.")
    cat("         They should  be of the dimension (n x n), and a vector of length n.")
    cat(" The dimensions are: \n")
    cat(" dim(P)      = ", dim(P), "\n")
    cat(" length(y)   = ", length(y), "\n")
    stop(call. = FALSE)  

  }



  ycolmat <- matrix(data=y, ncol=1)  ## makes it easier when dealing with this in Rcpp
  fnamebin <- fullpath("Mt.bin")
  if(!any(is.na(selected_loci))) selected_loci <- selected_loci-1
  ar <- calculate_reduced_a_rcpp(f_name_bin = fnamebin, varG=varG, P=P, 
                                 y=ycolmat, max_memory_in_Gbytes=availmemGb, 
                                 dims=dim_of_bin_M , selected_loci = selected_loci , 
                                 quiet = quiet )





## t(t(y)) is a trick to get y as a row matrix 
##return( varG * invMMt %*% P %*% t(t(y)))
return(ar)

}





calculate_a_and_vara <- function(maxmemGb=8, dims=NULL,
                         selectedloci = NA,
                         invMMtsqrt=NULL, transformed_a=NULL, transformed_vara=NULL,
                         quiet = 0,
                         indxNA = NULL)
{
 ## an Rcpp function to take dimension reduced a (BLUP) values 
 ## and transform them into the original a (BLUP) values and their variances 
 ## Args:
 ##   maxmemGb         maximum available memory (in Gigabytes) that are available for use
 ##   dims             a 2 element numeric vector with the number of rows,columns in M 
 ##   invMMtsqrt       a matrix object of the form (M %*% M^T)^{-0.5}
 ##   transformed_a    a numeric vector of the dimension reduced BLUP or a values
 ##   transformed_vara a numeric matrix of dimension dims(1) x dims(1) for the dimension reduced BLUPs (or a) values. 
 ##   selectedloci     an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                    have been selected to act as fixed QTL effects in the model. 



  file_bin <- fullpath("Mt.bin")
  if(!file.exists(file_bin)){
      cat("\n\n  Error: ", file_bin, " does not exist and it should have been created. \n\n")
      stop(call. = FALSE)
  }

  dimsMt <- c(dims[2], dims[1]) 

 
  if(!any(is.na(selectedloci))) selectedloci <- selectedloci-1
  calculate_a_and_vara_rcpp(f_name_bin=file_bin,
                    selected_loci = selectedloci,
                    inv_MMt_sqrt=invMMtsqrt,  
                    dim_reduced_vara = transformed_vara,
                    max_memory_in_Gbytes=maxmemGb, 
                    dims=dimsMt, 
                    a = transformed_a, 
                    quiet = quiet,
                    indxNA = indxNA) 

}


calculate_reduced_vara <- function(X=NULL, varE=NULL, varG=NULL, invMMt=NULL, MMtsqrt=NULL, quiet=FALSE)
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

  vars <- varG * diag(nrow(D1))  - ( D1 + D1 %*% C %*% solve(A - B %*% D1 %*% C) %*% B %*% D1 )

    return(vars )

}



check.inputs <- function(ncpu=NULL, availmemGb=NULL, 
                         file_genotype=NULL, 
                         file_phenotype=NULL )
{





if(!is.null(ncpu)){
 if(!is.numeric(ncpu)){
   cat("Error:  ncpu is not a numeric value. It is of class ", class(ncpu), "It should be the number of cpu.\n")
   return(TRUE)
 }
 if(ncpu < 1){
    cat("Error: ncpu cannot be a zero or a negative number. It should be the number of cpu. \n")
    return(TRUE)
 }
    
}

if(!is.null(availmemGb))
{
 if(!is.numeric(availmemGb)){
   cat("Error: availmemGb is not a numeric value. It is of class ", class(availmemGb), "It should be the number of giggabytes of RAM available. \n")
   return(TRUE)
 }
 if(availmemGb <= 0){
    cat("Error: availmemGb cannot be zero or a a negative number.  It should be the number of giggabytes of RAM available. \n")
    return(TRUE)
  } 
}


if(!is.null(file_genotype))
{
  genofile <- fullpath(file_genotype)

  if(!file.exists(genofile)){
    cat("Error: Cannot find marker file ", genofile, "\n")
    cat("       This could be a problem with the name of the file and/or the location of the file. \n")
    cat("       Perhaps specify the full name of the file (i.e. absolute directory path and file name) \n")
    cat("       Type help(ReadMarker) and go to the examples section for an example of this. \n") 
    return(TRUE)
  }
}


if(!is.null(file_phenotype))
{ 
  phenofile <- fullpath(file_phenotype)

  if(!file.exists(phenofile)){
    cat("Error: Cannot find phenotype file ", phenofile, "\n")
    cat("       This could be a problem with the name of the file and/or the location of the file. \n")
    cat("       Perhaps specify the full name of the file (i.e. absolute directory path and file name) \n")
    cat("       Type help(ReadPheno) and go to the examples section for an example of this. \n") 
    return(TRUE)
  }
}

  return(FALSE)


}

#' @title Read phenotype file
#' @description Read in the phenotypic data. The first row is assumed to contain the column headings. 
#' @param filename contains the name of the phenotype  file. The file name needs to be in quotes.
#' @param header a logical value. When \code{TRUE}, the first row of the file contains the names of the columns.  Default is \code{TRUE}.
#' @param csv   a logical value. When \code{TRUE}, a csv file format is assumed. When \code{FALSE}, a space separated format is assumed. Default
#'              is \code{FALSE}.
#' @param ... any other \code{\link{read.table}} formatting commands
#' @details  
#' 
#' \code{ReadPheno} reads in the phenotypic data. A space separated plain text file is assumed. Each row in this file 
#' corresponds to an individual. The number of rows in the phenotypic file must be the same as the number of rows in 
#' the marker data file. Also, the ordering of the individuals must be the same in the two files.  A space separated file with 
#' column headings is the default but can be changed with the \code{header} and \code{csv} options. 
#'
#' The phenotypic file may contain multiple traits and explanatory variables. 
#'
#' Missing values are allowed and should be coded NA. 
#
#' As an example, suppose we have three individuals for which we have collected data on two quantitative traits (y1 and y2), and 
#' four explanatory variables (var1, var2, var3, and var4). The data looks like  
#' \tabular{cccccc}{
#'     y1      \tab y2      \tab  var1  \tab var2 \tab   var3  \tab var4 \cr
#'     112.02  \tab -3.123  \tab  10    \tab H1   \tab   M     \tab blue \cr
#'     156.44  \tab 1.2     \tab  12    \tab H1   \tab   NA     \tab red  \cr
#'     10.3    \tab NA   \tab  9     \tab H2   \tab   F     \tab green 
#'}
#'
#' Then to load these data, we would use the command 
#'
#' \preformatted{pheno_obj <- ReadPhenr(filename="pheno.dat")}
#'
#' where "pheno.dat" is the name of the phenotype file, and \code{pheno_obj} is the R object that contains the 
#' results from reading in the phenotypic data.    
#' The column headings of the dataframe are "y1", "y2", 
#' "var1", "var2", "var3", and "var4" where "var2", "var3", and "var4" are read in as factors with the other columns treated as 
#' numerics. 
#'
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
#'
#' @seealso \code{\link{ReadMarker}} for reading in marker data, \code{\link{AM}} for performing association mapping.
#' @return 
#' a data frame is returned of the phenotypic data. If \code{header} is true, the 
#' names of the columns will be as specified by the first row of the phenotypic file. If \code{header} is \code{FALSE}, 
#' generic names are supplied by R in the form of V1, V2, etc.  If no column headings are given, these 
#' generic names will need to be used in the \code{trait} and \code{feffects} parameters in 
#' \code{\link{AM}}.  You can print out the column names of the dataframe by using
#'
#' \preformatted{names(pheno_obj)}
#'
#'
#'
#'
#'
#' @examples
#' # Read in  phenotypic data from ./extdata/
#' 
#' # find the full location of the phenotypic data 
#' complete.name <- system.file("extdata", "pheno.txt", package="AMplus")
#'
#' pheno_obj <- ReadPheno(filename=complete.name)
#'   
#'  ## print a couple of lines of the data file
#'  head(pheno_obj)
#'
ReadPheno <- function(filename = NULL, header=TRUE, csv=FALSE, ...){

  phenofile <- fullpath(filename)
  error.code <- check.inputs(file_phenotype=phenofile)
  if(error.code)
     stop(" ReadPheno has terminated with errors.", call. = FALSE)

  sep <- ""
  if(csv) sep=","
  phenos <- read.table(phenofile, header=header, sep=sep, ...)
cat("\n\n Loading Pheotype file ... \n\n")
cat("               Summary of Phenotype File  \n")
cat("              ~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
cat(" File name:                   ",  phenofile, "\n")
cat(" Number of individuals:       ", nrow(phenos), "\n")
cat(" Number of columns:           ", ncol(phenos), "\n\n")
cat(" First 5 rows of the phenotype file are \n")
print(head(phenos,n=5))
cat("\n Column classes are  \n")
for(ii in 1:ncol(phenos))
  cat(c( sprintf("%20s   %15s", names(phenos)[ii], class(phenos[[ii]]) ), "\n"))


cat("\n WARNING: if the column classes are incorrect, these will need to be changed by the user.\n\n\n")

  return(phenos)


}







#' @title Read map file
#' @description Read in the marker map  data. The first row is assumed to contain the column headings. 
#' @param filename contains the name of the map file. The file name needs to be in quotes.
#' @param csv   a logical value. When \code{TRUE}, a csv file format is assumed. When \code{FALSE}, a space separated format is assumed. 
#' @details
#' Unlike classical linkage mapping, association mapping does not require a map in order to analyse the data. 
#' So, reading in a map file is optional. If a map file is read into our package, then the marker names are used when 
#' results are reported from \code{\link{AM}}. If a map file is not supplied, then generic names M1, M2, ..., are given to the 
#' marker loci in columns 1, 2, ...
#' 
#' A space separated text file with column headings is assumed as the default input. The map file can have three or four 
#' columns. If the map file has three columns, then it is assmed that the three columns are the marker locus names, 
#' the chromosome number, and the map position (in any units). If the map file has four columns as with a PLINK map file, 
#' then the columns are assumed to be the marker locus names, the chromosome number, the map position in centimorgans, 
#' and the map position in base pairs. 
#'
#' Missing values are not allowed. 
#' 
#' The order of the marker loci in this file is assumed to be in the same order as the loci in the marker data file. 
#' @seealso \code{\link{ReadMarker}} and \code{\link{ReadPheno}}.
#' @return 
#' a data frame is returned of the map data. 
#'
#' @examples
#' # Read in  example map data from ./extdata/
#' 
#' # find the full location of the map data 
#' complete.name <- system.file("extdata", "map.txt", package="AMplus")
#'   
#' # read in map data 
#' map_obj <- ReadMap(filename=complete.name) 
#'                                
#'# look at first few rows of the map file
#' head(map_obj)
#'
#'
ReadMap  <- function( filename = NULL, csv=FALSE)
{
 mapfile <- fullpath(filename)
 error.code <-  check.inputs(file_genotype=filename)
 if(error.code)
    stop(" ReadMap has terminated with errors.", call. = FALSE)

  sep=""
  if(csv) sep=","
  map <- read.table(mapfile, header=TRUE, sep=sep)
cat("\n\n Loading map file ... \n\n")
cat("                    Summary of Map File  \n")
cat("                   ~~~~~~~~~~~~~~~~~~~~~~ \n")
cat(" File name:                   ",  mapfile, "\n")
cat(" Number of marker loci:       ", nrow(map), "\n")
cat(" Number of columns:           ", ncol(map), "\n")
cat(" Number of chromosomes:       ", length(unique(map[[2]])), "\n\n")
cat(" First 5 markers of the map file are \n")
print(head(map, n=5))
cat("\n\n")

return(map)

}



create.bin  <- function(file_genotype=NULL,  type="text", AA=NULL, AB=NULL, BB=NULL, 
                         availmemGb=8, dim_of_bin_M=NULL, csv=FALSE, quiet=FALSE){
 ## an Rcpp function to create the packed binary file of the genotype data M and Mt
 ## from marker data. The marker data may be from an ASCII file or PLINK ped file.
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## AA, AB, BB       numeric codes for associated genotypes in marker genotype file
 ## availmemGb     available memory for converstion to packed binary
 ## dim_of_bin_M             row, column dimensions of M.  
 ## type            where file type is text or PLINK
 binMfile <- fullpath("M.bin")
 binMtfile <- fullpath("Mt.bin")


if (type=="text"){
    ## text genotype file
    createM_rcpp(f_name = file_genotype, type=type ,  f_name_bin = binMfile, AA = AA, AB = AB, BB = BB,
               max_memory_in_Gbytes=availmemGb,  dims = dim_of_bin_M , csv = csv,
               quiet = quiet)
    createMt_rcpp(f_name = file_genotype, f_name_bin = binMtfile,  AA = AA, AB = AB, BB = BB,
                  max_memory_in_Gbytes=availmemGb,  dims = dim_of_bin_M, csv=csv, quiet = quiet )
} else {
    ## PLINK ped file
    ## using -9 to indicate missing/null genotypes
    createM_rcpp(f_name = file_genotype, type=type,  f_name_bin = binMfile, AA ="-9", AB = "-9", BB = "-9",
               max_memory_in_Gbytes=availmemGb,  dims = dim_of_bin_M , csv=csv, quiet = quiet)
    createMt_PLINK_rcpp(f_name = file_genotype, f_name_bin = binMtfile,   
                  max_memory_in_Gbytes=availmemGb,  dims = dim_of_bin_M, quiet = quiet )


}  ## end if else type

 return(NULL)

}




#' @title Read  marker data.
#' 
#' @description
#' A function for reading in marker data. Three different types of data can be read. 
#' @param filename contains the name of the marker  file. The file name needs to be in quotes. 
#' @param type  specify the type of file. Choices are "text" and "PLINK".
#' @param AA     the character(s) or number corresponding to the AA genotype in the marker genotype file. 
#' This must be specified if the file type is "text".  The character(s) must be in quotes.
#' @param AB     the  character(s) or number  corresponding to the AB genotype in the marker genotype file. This can be left unspecified 
#'               if there are no heterozygote genotypes (i.e. the individuals are inbred). If specified, the character(s) must 
#'               be in quotes. 
#' @param BB        the character(s) or number corresponding to the BB genotype in the marker genotype file. This must be specified if the file type is "text".  The character(s) must be in quotes.
#' @param availmemGb a numeric value. It specifies the amount of available memory (in Gigabytes). 
#'         This should be set to be as large as possible.  
#' @param csv   a logical value. When \code{TRUE}, a comma separated text file is assumed. When \code{FALSE}, a space separated text file 
#'         is assumed. Only type="text" files are permitted to be space or comma separated.  
#'         PLINK ped files are assumed to be space separated.
#' @param  quiet      an integer value specifying the number of marker loci for which diagnostic information is 
#' to be printed to the screen. This is useful for error checking. 
#'
#' @details
#' 
#' \code{ReadMarker} can handle three different types of marker data; namely
#' previously read marker data, genotype data in a plain text file, and PLINK ped files. 
#'  
#' \subsection{\strong{Reading in previously read marker data}}{
#' To read marker data that has been previously read with \code{ReadMarker} in another R session, run 
#' the function with no arguments. For example 
#'
#' \preformatted{geno_obj <- ReadMarker()}
#'
#' where \code{geno_obj} is the name of the user defined R object that is to contain the marker data. This can be 
#' anything you want. 
#'
#' For this command to work without error, the working directory needs to 
#' be the same as the working directory form which the \code{\link{ReadMarker}} was first run. 
#' You can check on what the 
#' current working directory is with the command 
#' \preformatted{getwd()}
#' and you can change your working directory with the command 
#' \preformatted{changewd("\\newworkingdir")}
#' where "\\newworkingdir" is the new working directory for the R session.
#'
#' Loading the marker data in this way is much faster than reading the data again from file 
#' because there is no need to pre-process the data, nor check the data for errors. }
#'
#' \subsection{\strong{Reading in plain text file containing marker genotypes}}{
#' To load a text file that contains snp genotypes, run \code{ReadMarker} with \code{filename} set to the name of the file, 
#' and \code{AA}, \code{AB}, \code{BB} set to corresponding genotype values. 
#' The genotype values in the text file can be numeric, character, or even a mix of both.
#'
#' We make the following assumptions
#' \itemize{
#' \item{The text file does not contain row or column headings}
#' \item{Individuals are diploid}
#' \item{The rows of the text file are the individuals and the columns are the marker loci}
#' \item{The file is either space or comma separated}
#' \item{The mapping of the observed genotypes to \code{AA}, \code{AB}, and \code{BB}, remains the same for all loci}
#' \item{Individuals are outbred when \code{AA}, \code{AB}, and \code{BB} are specified and 
#'  inbred when only \code{AA}, and \code{BB} are specified}
#'}
#'
#' Suppose we have a space separated text file with marker genotype data collected from five snp loci on three individuals
#' where the snp genotype AA has been coded 0, the snp genotype AB has been coded 1, and the snp genotype BB has been coded 2
#' \tabular{ccccc}{
#'  0 \tab  1  \tab  2 \tab  0\tab   2 \cr
#'  1 \tab  1  \tab  0 \tab  2 \tab  0 \cr
#'  2 \tab  2  \tab  1 \tab  1 \tab  1
#'}
#'  Then to load these data, we would use the command
#'
#' \preformatted{geno_obj <- ReadMarker(filename="geno.txt", AA=0, AB=1, BB=2)}
#'
#' where the results from running the function are placed in \code{geno_obj}, and the name of the file is "geno.txt" 
#' which is located in the current working directory (which can be found 
#' by using the command \code{getwd()}).
#'
#' Suppose we have a comma separated text file with marker genotype data collected from five snp loci on three individuals 
#' where the snp genotype AA has been coded a/a, the snp genotype AB has been coded a/b, and the snp genotype BB has been coded b/b
#' \tabular{c}{
#'  a/a,a/b,b/b,a/a,b/b \cr
#'  a/b,a/b,a/a,b/b,a/a \cr
#'  b/b,b/b,a/b,a/b,a/b
#'}
#'
#' Then to load these data, we would use the command 
#'
#' \preformatted{geno_obj <- ReadMarker(filename="geno.csv", AA="a/a", AB="a/b", BB="b/b", csv=TRUE)}
#'
#' where \code{geno_obj} is used by \code{\link{AM}}. 
#'}
#'
#' \subsection{\strong{Reading in a PLINK ped file}}{
#' PLINK is a well known toolkit for the analysis of genome-wide association data. See  \url{http://pngu.mgh.harvard.edu/~purcell/plink/}
#' for details. 
#'
#' Full details of PLINK ped files can be found \url{http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped}. Briefly, 
#' the PED file is a space delimited file (tabs are not allowed): the first six columns are mandatory:
#'
#' \tabular{l}{
#'     Family ID  \cr
#'     Individual ID \cr
#'     Paternal ID \cr
#'     Maternal ID \cr
#'     Sex (1=male; 2=female; other=unknown) \cr
#'     Phenotype
#'}
#'
#'
#' Here, these columns can be any values since \code{ReadMarker} ignores these columns.  
#'
#' Genotypes (column 7 onwards) can be any character 
#' (e.g. 1,2,3,4 or A,C,G,T or anything else) except 0 which is, by default, the missing genotype character. All markers should be biallelic. 
#' All SNPs must have two alleles specified.  Missing alleles (i.e 0 or -) are not allowed and will generate an error. 
#' No column headings should be given. 
#' 
#' As an example, suppose we have data on three individuals  genotyped for four snp loci 
#' \tabular{cccccccccccccc}{
#'     FAM001 \tab 101  \tab  0    \tab 0   \tab   1  \tab 0 \tab  A  \tab  G  \tab C \tab C \tab C \tab G \tab A \tab A \cr
#'     FAM001 \tab 201  \tab  0    \tab 0   \tab   2  \tab 0 \tab  A  \tab  A  \tab C \tab T \tab G \tab G \tab T \tab A \cr 
#'     FAM001 \tab 300  \tab  101  \tab 201 \tab   2  \tab 0 \tab  G  \tab  A  \tab T \tab T \tab C \tab G \tab A \tab T 
#'}
#'
#' Then to load these data, we would use the command 
#'
#' \preformatted{geno_obj <- ReadMarker(filename="PLINKfile.ped", type="PLINK")}
#'
#' where \code{geno_obj} is used by \code{\link{AM}}. 
#'}
#'
#'
#' \subsection{Reading in other formats}{
#' It is possible to convert other file formats into PLINK ped files. See \url{http://pngu.mgh.harvard.edu/~purcell/plink/} for details. 
#' For example, to convert  vcf file into a PLINK ped file, use the PLINK command
#' 
#' \preformatted{PLINK --vcf filename.vcf --recode --out newfilename}
#'
#' and to convert a binary ped file (bed) into a ped file, use the PLINK command
#'
#' \preformatted{PLINK --bfile filename --recode --tab --out newfilename}
#'  
#'}
#'
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
#' proportion of missing data is low.  See  George and Cavanagh (2015) for details.  
#' }
#'
#' @return  To allow \code{\link{AM}} to handle data larger than the memory capacity of a machine, \code{ReadMarker} doesn't load 
#' the marker data into memory. Instead, it creates a packed binary file of the marker data and its transpose. The object returned by
#' \code{ReadMarker} is a list object with the elements \code{binfileM} , \code{binfileMt}, and \code{dim_of_bin_M}  
#' which is the full file name (name and path)  of the packed binary file for the marker  data,  the full file name of the packed binary file 
#' for the transpose of the marker  data,  and a 2 element vector with the first element the number of individuals and the second 
#' element the number of marker loci. 
#' 
#' @references George AW and Cavanagh C. 2015. Genome-wide Association Mapping in Plants. 
#' Theorectical and Applied Genetics 128: 1163-1174.
#'
#' @examples
#'   #--------------------------------
#'   #  Example 1
#'   #-------------------------------
#'   #
#'   # Read in the genotype data contained in the text file "geno.txt"
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file("extdata", "geno.txt", package="AMplus")
#'
#'   
#'   # Here,  0 values are being treated as genotype AA,
#'   # 1 values are being treated as genoytpe AB, 
#'   # and 2 values are being treated as genoytpe BB. 
#'   # 4 Gbytes of memory has been specified. The file is space separated with the rows the individuals
#'   # and the columns the snp loci.
#'   geno_obj <- ReadMarker(filename=complete.name, type="text", AA=0, AB=1, BB=2, availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
#'   #--------------------------------
#'   #  Example 2
#'   #-------------------------------
#'   #
#'   # Read in the allelic  data contained in the PLINK ped file "geno.ped"
#'   #
#'   # The function system.file() gives the full file name (name + full path).
#'   complete.name <- system.file("extdata", "geno.ped", package="AMplus")
#'
#'   
#'   # Here,  the first 6 columns are being ignored and the allelic 
#'   # information in columns 7 -  10002 is being converted into a packed binary file. 
#'   # 4 Gbytes of memory has been specified. The file is space separated with the rows the individuals
#'   # and the columns the snp loci.
#'   geno_obj <- ReadMarker(filename=complete.name, type="PLINK", availmemGb=4) 
#'    
#'   # view list contents of geno_obj
#'   print(geno_obj)
#'
#'
ReadMarker <- function( filename=NULL, type="text",
                           AA=NULL, AB=NULL, BB=NULL, 
                           availmemGb=8, 
                           csv = FALSE,
                           quiet=FALSE){


 if (nargs() == 0){
   ## if no arguments are supplied to ReadMarker, then it is assumed that ReadMarker has been run
   ## previously and the file M.RData  is in the current working directory. Need to check this. 
   ## checks if M.RData which contains list object geno, Mt.bin, and M.bin exist. IF so, it returns list object geno

   if(file.exists(fullpath("M.RData")))
   {
       # check that M.bin and Mt.bin exist in this directory
       if(!file.exists(fullpath("M.bin"))){
         cat(" The binary file M.bin could not be found in current working directory ", getwd(), "\n")
         cat(" This file is created when ReadMarker is run with either a text file or PLINK ped file as input. \n")
         cat(" Supply a file name to ReadMarker. Type  help(ReadMarker) for more detals \n")
         stop(" ReadMarker has terminated with errors \n", call. = FALSE)
       }

       if(!file.exists(fullpath("Mt.bin"))){
         cat(" The binary file Mt.bin could not be found in current working directory ", getwd(), "\n")
         cat(" This file is created when ReadMarker is run with either a text file or PLINK ped file as input. \n")
         cat(" Supply a file name to ReadMarker. Type  help(ReadMarker) for more detals \n")
         stop(" ReadMarker has terminated with errors \n", call. = FALSE)
       }
       ## looks like everything is good. Return geno list object
       cat(" The files M.RData, M.bin, and Mt.bin, in current working directory ", getwd(), " have been found and will be used for the association mapping analysis. \n")
   load("M.RData")
       return(geno)

   } else {

       cat(" The R file M.RData could not be found in current working directory ", getwd(), "\n")
       cat(" This file is created when ReadMarker is run with either a text file or PLINK ped file as input. \n")
       cat(" Supply a file name to ReadMarker. Type  help(ReadMarker) for more detals \n")
       stop(" ReadMarker has terminated with errors \n", call. = FALSE)

   }  ## end if(file.exists(fullpath("M.RData"))

 }  else {
      ## read in either a text file or a PLINK file. The paramter type must be specified. Default it text file. 

   if(is.null(type)){
      cat(" type must be set to \"text\" or \"PLINK\". \n")
      stop(" ReadMarker has terminated with errors \n", call. = FALSE)
   }
   if(!(type=="text" || type=="PLINK") ){
      cat(" type must be set to \"text\" or \"PLINK\". \n")
      stop(" ReadMarker has terminated with errors \n", call. = FALSE)
   }

    ## ------   PLINK ped file -------------------
    if (type=="PLINK"){

       ## checking if a PLINK file has been specified. 
       if (is.null(filename)){
            cat(" The name of the PLINK ped file is missing in current working directory ", getwd(), " \n")
            cat("  Run ReadMarker with filename  set to the name of the PLINK ped input file and set type=\"PLINK\". \n")
            stop(" ReadMarker has terminated with errors \n", call. = FALSE)
       }
       if (!file.exists(fullpath(filename) )){
            cat(" The PLINK ped file ", filename, " could not be found in the current directory ", getwd(), "\n")
            stop(" ReadMarker has terminated with errors \n", call. = FALSE)
       }

       ## Rcpp function to get dimensions of PLINK ped  file
       dim_of_bin_M <- getRowColumn(fname=fullpath(filename), csv=FALSE )

       ## Rcpp function to create binary packed M and Mt file 
       create.bin(file_genotype=fullpath(filename), type=type, availmemGb=availmemGb, dim_of_bin_M=dim_of_bin_M,  quiet=quiet  )
       binfileM <- fullpath("M.bin")
       binfileMt <- fullpath("Mt.bin")


   }  else {
      ## ------------  text file -----------------------
      ## Assuming a text file that may be comma separated with numeric genotypes that need to be mapped onto AA, AB, and BB. 
      ## check of parameters
      error.code <- check.inputs(file_genotype=filename, availmemGb=availmemGb)
      if(error.code)
          stop(" ReadMarker has terminated with errors.", call. = FALSE)
 ## Has AA, AB, BB been assigned character values
  if(is.null(AA) ||  is.null(BB))
  { 
     cat("Error: The function parameters AA and BB must be assigned a numeric or character value since a text file is being assumed. \n")
     cat("       Type help(ReadMarker) for help on how to read in the marker data. \n")
     stop(" ReadMarker has terminated with errors \n", call. = FALSE)
  }

  ## if there are no hets. 
  if(is.null(AB))
     AB <- "NA"  ## no hets 




  genofile <- fullpath( filename) 



  ## Rcpp function to get dimensions of ASCII genotype file
  dim_of_bin_M <- getRowColumn(fname=genofile, csv=csv )

  ## Rcpp function to create binary packed M and Mt file from 
  create.bin(file_genotype=genofile, type=type, AA=as.character(AA), AB=as.character(AB), BB=as.character(BB), 
              availmemGb=availmemGb, dim_of_bin_M=dim_of_bin_M, csv=csv, quiet=quiet  )
    binfileM <- fullpath("M.bin")
    binfileMt <- fullpath("Mt.bin")

 
  }  ## end if else nargs()==1  (PLINK case)

  if(type=="PLINK"){
     ## PLINK ped file which has been converted into a geno file
     ## need to adjust dim_of_bin_M[2]  -- columns
     dim_of_bin_M[2]<- (dim_of_bin_M[2] - 6)/2  ## allelic info starts in col 8
  }

  geno <- list("binfileM"=binfileM, "binfileMt"=binfileMt,
               "dim_of_bin_M" = dim_of_bin_M)
  save(geno, file="M.RData")

  ## create M.Rdata file in current directory
  return(geno)

  } ## end if else nargs()==0

}  ## end function call ReadMarker




extract_geno <- function(colnum=NULL, availmemGb=8, 
                          dim_of_bin_M=NULL, 
                          selected_locus=NA,
                          indxNA = NULL )
  {
    ## Rcpp function to extra a column of genotypes from a  binary packed file of M

    binfileM <- fullpath("M.bin")
    selected_locus <- colnum - 1  ## to be consistent with C++'s indexing starting from 0
    ## AWG if(!any(is.na(indxNA))) indxNA <- indxNA - 1
    if (!(length(indxNA)==0)) indxNA <- indxNA - 1

    geno <- extract_geno_rcpp(f_name_bin=binfileM, max_memory_in_Gbytes = availmemGb, 
                              selected_locus=selected_locus, dims=dim_of_bin_M,
                              indxNA = indxNA)

    return(geno)

  }



constructX <- function(currentX=NULL, loci_indx=NULL, 
                       availmemGb=8, dim_of_bin_M=NULL,
                       indxNA = NULL , map=NULL)
  {
    ## R function to construct the design matrix X
    ## Args
    ##   currentX    current model matrix
    ##   loci        the marker loci to be included as fixed QTL effects (additive model)
    ##   indxNA      those individuals that should be removed due to missing phenotypes
   
   if(is.na(loci_indx))
   {
     return(currentX)
   } else {
       genodat <- extract_geno(colnum=loci_indx, 
                           availmemGb=availmemGb, dim_of_bin_M=dim_of_bin_M,
                           indxNA = indxNA)
      newX <- cbind(currentX, genodat)
      colnames(newX) <- c(colnames(currentX), as.character(map[[1]][loci_indx])) ## adding col names to new X  
      return(newX)
   }
  }














