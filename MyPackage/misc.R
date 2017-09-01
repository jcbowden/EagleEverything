

ReshapeM  <- function(fnameM, fnameMt, indxNA, dims){
   ## function to create a temp version of M.ascii and Mt.ascii where the rows and columns, 
   ## respectively have been removed for the elements in indxNA
   
   ## its indxNA-1 so that indexes start from 0 as in c++
   res <- ReshapeM_rcpp(fnameM=fnameM, fnameMt=fnameMt, indxNA=(indxNA-1), dims=dims)
   return(res)  ## returns integer vector with new dims of reshaped matrix M
}


doquiet <- function(dat, num_markers, lab){
     ## a diagnostic function for printing the contents of matrix or vector
     ## used for error checking

     if(dim(dat)[1] == 1 || dim(dat)[2] ==1 )
         dat <- as.numeric(dat)

     if(class(dat)=="matrix"){
          ### if dat is a matrix

         if(num_markers > 0){
           message(" Dimension of ", lab, " is ", dim(dat)[1], " by ", dim(dat)[2], " \n")
           message(" First few rows and ", num_markers, " columns of ", lab, " are: \n")
           if(nrow(dat) > 5 && ncol(dat) > num_markers){
             for(xx in 1:5)
               message(sprintf(" %f ", dat[xx, 1:num_markers]))

           }
           if(nrow(dat) <=5  &&  ncol(dat) > num_markers)
             for(xx in 1:nrow(dat))
               message(sprintf(" %f ", dat[xx, 1:num_markers]))
           if(nrow(dat) > 5  &&  ncol(dat) <=  num_markers)
             for(xx in 1:5)
               message(sprintf(" %f ", dat[xx, 1:ncol(dat)]))
           if(nrow(dat) <= 5  &&  ncol(dat) <=  num_markers)
             for(xx in 1:nrow(dat))
               message(sprintf(" %f ", dat[xx, 1:ncol(dat)]))
           message("\n\n")
         }
     } ## end if class(dat)

     if(class(dat)=="numeric" || class(dat)=="vector"){
       if(num_markers > 0){
          message(" Length of ", lab, " is ", length(dat), "\n")
          message(" The first ", num_markers, " elements of the vector are ", lab, "\n")
          if(length(dat) > num_markers)
             message(sprintf(" %f ", dat[1:num_markers]))
          if(length(dat) <= num_markers)
             message(sprintf(" %f ", dat[1:length(dat)]))
       message("\n\n")
       }
    }

    if(!(class(dat)=="matrix" || class(dat)=="vector" || class(dat)=="numeric"))
      message(" Internal error in function doquiet. dat not matrix or vector or numeric. \n")

}



.form_results <- function(trait, selected_loci, map,  fformula, indxNA,
                           ncpu, availmemGb, quiet,  extBIC )
{
  ## internal function - used by AM for forming the results object
  if (length(selected_loci) > 1){
   sigres <- list(trait=trait,
                    fformula = fformula,
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
                    fformula = fformula,
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
    ## internal function: use only in AM function
    ## title
    message("\n\n\n")
message("                    Multiple-Locus Association Mapping")
message("                            Version 1.0 \n")
message(" ")
message("   . ,-\"-.   ,-\"-. ,-\"-.   ,-\"-. ,-\"-. ,-\"-. ,-\"-.   ,-\"-. ,-\"-.    ")  
message("    X | | \\ / | | X | | \\ / | | X | | \\ / | | X | | \\ / | | X | | \\ /   ")
message("   / \\| | |X| | |/ \\| | |X| | |/ \\| | |X| | |/ \\| | |X| | |/ \\| | |X|   ")
message("      `-!-' `-!-\"   `-!-' `-!-'   `-!-' `-!-\"   `-!-' `-!-'   `-!-' `-     \n\n")

}




.build_design_matrix <- function(pheno=NULL,  indxNA=NULL, fformula=NULL, quiet=TRUE  )
{
   ## internal function: use only in AM function and SummaryAM  function
   ## build design matrix given character vector fformula of column names

   ## assign model matrix X
   if(is.null(fformula))
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
        Xmat <- model.matrix(fformula, data=pheno)
     }  else {
        # there is an issue with creating Xmat when it includes
        # factors that have some of their levels removed. 
        ph <- pheno[-indxNA,]
        mat <- get_all_vars(formula=fformula, data=ph)
        for(ii in names(mat)){
           if(is.factor(ph[,ii])){
              ph[,ii] <- as.factor(as.character(ph[,ii]))
           }
        }  ## for    
        Xmat <- model.matrix(fformula, data=ph)
     } ## if else (length(indxNA)==0)
   } 

 if (!quiet ){
   message("Dimension of design matrix, before addition of marker fixed effects is ", nrow(Xmat), "rows and ", ncol(Xmat), "columns.\n") 
 }

if(!is.matrix(Xmat))
   Xmat <- matrix(data=Xmat, ncol=1)

## remove column that are 0 sum 
indx <- which(colSums(Xmat)==0)
if(length(indx) > 0)
   Xmat <- Xmat[, -indx]


  return(Xmat)
}


.calcMMt <- function(geno, availmemGb, ncpu, selected_loci, quiet)
  {
    ## internal function: used only in AM  and SummaryAM
    ## calculates M %*% t(M) via C++ for out of memory calculation
    MMt <- calculateMMt(geno=geno[["asciifileM"]], availmemGb=availmemGb, 
                           ncpu=ncpu, 
                           dim_of_ascii_M = geno[["dim_of_ascii_M"]], 
                           selected_loci=selected_loci, quiet = quiet, message=message) 
    gc()


    ## Trick for dealing with singular MMt due to collinearity
    MMt <- MMt/max(MMt) + diag(0.95, nrow(MMt)) 
    return(MMt)
  }


  .calcVC <- function(trait, currentX, MMt, ngpu)
  {
    ## internal function: used by AM and SummaryAM
    ## perform likelihood ratio test for variance component Var_g
    res_full <- emma.REMLE(y=trait, X= currentX , K=MMt, ngpu=ngpu)
    return(list("vg"=res_full$vg, "ve"=res_full$ve))

  }

 .calc_extBIC <- function(trait=NULL, currentX=NULL, MMt=NULL,  geno=NULL, quiet=TRUE)
 { 
   ## internal function: used by AM 
   ## smallest extBIC and BIC is best
   ## internal function: use in AM only
   res_p <- emma.MLE(y=trait, X= currentX , K=MMt, llim=-100,ulim=100)
   BIC <- -2 * res_p$ML + (ncol(currentX)+1) * log(length(trait))  ## fixed effects + variance component

   extBIC <- BIC + 2 * lchoose(geno$dim_of_ascii_M[2], ncol(currentX) - 1)  

    return(extBIC)
 }



 .print_header <- function(){
   ## internal function: used by AM
   message("\n\n\n                           Final  Results  \n")
   message(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
 }

.print_final  <- function(selected_loci, map,  extBIC )
{
  ## internal function: used by AM
  if (length(selected_loci) == 1 & any(is.na(selected_loci)))
  {
      message("No significant marker-trait associations have been found. \n\n")
  }  else {
     .print_results(selected_loci=selected_loci, map=map,  extBIC=extBIC)
          message("\n\n")
  }   ## end if else


}  ## end function print.finals

 .print_results <- function(itnum=NULL, selected_loci, map, extBIC)
 ## internal function: used by AM
 {  if(!is.null(itnum)){ 
       message(" Significant marker-trait association found. \n")
       message(" New results after iteration ", itnum, "are \n")
    }
    message(sprintf("%15s  %10s        %10s     %10s        %10s ", 
                 "SNP", "Chrm", "Map Pos",  "Col Number",       "extBIC"))
    message(sprintf("%15s  %10s        %10s     %10s        %10s ", 
                 "-----", "------", "---------",  "-----------",       "---------"))

    for(ii in 1:length(selected_loci)){
       if(is.na(selected_loci[ii])){
       message(sprintf("%15s  %10s        %10s        %8s           %-8.2f ", 
        "Null Model", " ", " ", " ", extBIC[ii] ))
       }  else {
       message(sprintf("%15s  %10s        %10s       %8s            %-8.2f ", 
        map[[1]][selected_loci[ii]], map[[2]][selected_loci[ii]], as.character(map[[3]][selected_loci[ii]]), 
             selected_loci[ii], extBIC[ii] ))
     }  ## end if else 
   }
    message("\n\n\n")
 }



  .find_qtl <- function(geno, availmemGb,  selected_loci, MMt, invMMt, best_ve, best_vg, 
                       currentX,  ncpu, quiet, trait, ngpu )
  {
    ##  internal function: use by   AM
    H <- calculateH(MMt=MMt, varE=best_ve, varG=best_vg, message=message ) 
    if(!quiet)
        doquiet(dat=H, num_markers=5, lab="H")

    P <- calculateP(H=H, X=currentX , message=message) 
    if(!quiet)
        doquiet(dat=P, num_markers=5 , lab="P")
    rm(H)
    gc()
 
    
    ## Looks at the stability of the MMt calculation especially if there are near identical rows of data in M
    error_checking <- FALSE
    if (!quiet )
       error_checking <- TRUE
    MMt_sqrt_and_sqrtinv  <- calculateMMt_sqrt_and_sqrtinv(MMt=MMt, checkres=error_checking, 
                              ngpu=ngpu , message=message) 
    if(!quiet){
       doquiet(dat=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], num_markers=5, lab="sqrt(M %*% M^t)")
       doquiet(dat=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]], num_markers=5, lab="sqrt(M %*% M^t)^-1")
    }
    if(!quiet ){
      message(" quiet =", quiet, ": beginning calculation of the BLUP estimates for dimension reduced model. \n")
    }
    hat_a <- calculate_reduced_a(varG=best_vg, P=P, 
                       MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                       y=trait, quiet = quiet , message=message)   
    if(!quiet)
       doquiet(dat=hat_a, num_markers=5, lab="BLUPs")


     rm(P)
     gc()

    if(!quiet ){
      message(" quiet = ", quiet, ": beginning calculation of the standard errors  of BLUP estimates for dimension reduced model. \n")
    }

    var_hat_a    <- calculate_reduced_vara(X=currentX, varE=best_ve, varG=best_vg, invMMt=invMMt, 
                                                MMtsqrt=MMt_sqrt_and_sqrtinv[["sqrt_MMt"]], 
                                                quiet = quiet, message=message ) 
    if(!quiet)
             doquiet(dat=var_hat_a, num_markers=5, lab="SE of BLUPs")


   
     gc()
    if(!quiet ){
      message(" quiet = ", quiet, ": beginning calculation of BLUPS and their standard errors for full model. \n")
    }

     a_and_vara  <- calculate_a_and_vara(geno = geno,
                                         maxmemGb=availmemGb, 
                                            selectedloci = selected_loci,
                                            invMMtsqrt=MMt_sqrt_and_sqrtinv[["inverse_sqrt_MMt"]],
                                            transformed_a=hat_a, 
                                            transformed_vara=var_hat_a,
                                            quiet=quiet, message=message) 
     if(!quiet){
        doquiet(dat=a_and_vara[["a"]], num_markers=5, lab="BLUPs for full model")
        doquiet(dat=a_and_vara[["vara"]], num_markers=5, lab="SE of BLUPs for full model")
     }

  
    ## outlier test statistic
    if (!quiet ) 
        message(" quiet = ", quiet, ": beginning calculation of outlier test statistics. \n")
    tsq <- a_and_vara[["a"]]**2/a_and_vara[["vara"]]
    if(!quiet)
       doquiet(dat=tsq, num_markers=5, lab="outlier test statistic")


    indx <- which(tsq == max(tsq, na.rm=TRUE))   ## index of largest test statistic. However, need to account for other loci 
                                         ## already having been removed from M which affects the indexing

    ## taking first found qtl
    indx <- indx[1]

    orig_indx <- seq(1, geno[["dim_of_ascii_M"]][2])  ## 1:ncols
    return(orig_indx[indx])
}



## internal function to get absolute path name under Unix and windows
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


##---------------------------------------------------------
##  The following emma functions are part of the EMMA package
##-----------------------------------------------------------
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
#    if(ngpu > 0){
#      if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
#         eig <- rcppMagmaSYEVD::eigen_mgpu(K, symmetric=TRUE)
#       }

#     } else {
      eig <- eigen(K, symmetric = TRUE)
#     }
    return(list(values = eig$values, vectors = eig$vectors))
}

emma.eigen.R.wo.Z <-  function (K, X, ngpu=0) 
{
    n <- nrow(X)
    q <- ncol(X)
    dn <- diag(n)
    S <- dn - X %*% solve(crossprod(X, X)) %*% t(X)
    gc()
#    if(ngpu > 0){
#     if(requireNamespace("rcppMagmaSYEVD", quietly = TRUE)) {
#       eig <- rcppMagmaSYEVD::eigen_mgpu(S %*% (K + dn) %*% S, symmetric = TRUE, only_values=FALSE)
#     }
#    } else {
       eig <- eigen(S %*% (K + dn) %*% S, symmetric = TRUE)
#    }


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


check.for.NA.in.trait <- function(trait=NULL)
{
     ## internal function for AM 
     ## to return the positions of NA in a trait
     ## ordered for largest to smallest (this is important for ReshapeM_rcpp code

       ## check for NA's in trait
        indxNA <- which(is.na(trait))
        if(length(indxNA)==0){
          indxNA <- vector("numeric", 0)
        } else {
          ## place in reverse order
          indxNA <- sort(indxNA, decreasing = TRUE)
message(cat("\n\n WARNING!!!! The individuals in rows ", indxNA, " either have missing trait data "))
message("             and/or missing explanatory variable values. These individuals have ")
message(cat("             been removed from the analysis.  \n"))
          if(any(is.na(indxNA))){
            message("Error:  (internal).  indxNA contains NA values. ")
            message(" AM has terminated with errors. ")
            return(NULL)
          }
        }

      return(indxNA)
   } 


check.inputs.mlam <- function (ncpu, availmemGb, colname.trait, map, pheno, 
                  geno )
{
  ## internal function for AM


if(is.null(colname.trait)){
   message("Error: the name of the column containing the trait data must be given. \n")
   return(TRUE)
}

if(is.null(pheno)){
   message("Error: the pheno parameter has not been set. ")
   message("       Set this parameter to the object that contains ")
   message("       the phenotype data. This object is the result of running  ")
   message("       ReadPheno. \n")
   return(TRUE)
}

if(is.null(geno)){
   message("Error: the marker data has not been specified.")
   message("       If you are using the GUI, then go back and read in your marker data. ")
   message("       If you are using 'Eagle' functions, then set the geno parameter to the ")
   message("       output from running ReadMarker(). ")
   return(TRUE)
}


if(class(try(class(geno), silent=TRUE)) == "try-error"){
   message("Error: the object supplied to the geno parameter does not exist. ")
   message("       This object is set by running ReadMarker. Type help(ReadMarker) for help ")
   message("       on running this command. \n")
   return(TRUE)
}


if(class(try(class(pheno), silent=TRUE)) == "try-error"){
   message("Error: the object supplied to the pheno parameter does not exist. ")
   message("       This object is set by running ReadPheno. Type help(ReadPheno) for help. ")
   message("       on running this command. \n")
   return(TRUE)
}


## checking list structure of geno
if(!is.list(geno)){
  message("Error: the geno object is not a list object. ")
  message("     The geno object is obtained from running ReadMarker.Type help(ReadMarker) for help. \n")
  return(TRUE)
}

## checking if pheno is a data frame 
if(!is.data.frame(pheno)){
  message("Error: the pheno object is not a data frame. ")
  message("      It is a ", class(pheno), "\n")
  return(TRUE)
}



nms <- names(geno)
indx <- match(nms, c("asciifileM", "asciifileMt", "dim_of_ascii_M" ))
if(any(is.na(indx))){
  message("Error: there is a problem with the list structure of the geno object. ")
  message("       It should contain the elements asciifileM, asciifileMt, and dim_of_ascii_M. ")
  message("       The object supplied contains the elements ", names(geno) )
  return(TRUE)
}

if(is.null(map)){
    message("WARNING: no map object has been specified. A generic map ")
    message("         will be assumed.                                ")
    map <- data.frame(Mrk= paste("M", 1:geno[["dim_of_ascii_M"]][2]), Chrm=1, Pos=1:geno[["dim_of_ascii_M"]][2])
}






 ## checks for colname.trait
 if(is.null(colname.trait)){
    message("Error: the column name for the trait/response has not been specified.")
    message("       Please set trait to the column name of the trait data in ")
    message("       the phenotype file. The allowable column names are ", names(pheno) )
    return(TRUE)
 }

 if(length(colname.trait)>1){
    message("Error: multiple column names for the trait have been specified. ")
    message("       Only a single column name should be  assigned to trait. ")
    return(TRUE)
 }

 indx <- match(colname.trait, names(pheno))
 if(any(is.na(indx))){
   message("Error: the trait column name does not match any of the column names in the phenotype file. ")
   message("       The name that has been supplied is ", colname.trait)
   message("       The column names of the phenotype file are ", names(pheno))
   return(TRUE)
 }







 ## check that geno and pheno contain the same number of individuals
 if(nrow(pheno) !=  geno[["dim_of_ascii_M"]][1])
 {
   message("Error: the number of individuals specified in the phenotype file is ", nrow(pheno))
   message("       the number of individuals specified in the genotype file is ",  geno[["dim_of_ascii_M"]][1])
   message("       The number of individuals should be the same in the two files.")
   return(TRUE)
 }

 ## check that map and geno contain the same number of snp
 if(nrow(map) != geno[["dim_of_ascii_M"]][2])
 {
   message("Error: the number of marker loci in the map file is ", nrow(map))
   message("       The number of marker loci in the genotype file is ", geno[["dim_of_ascii_M"]][2])
   message("       The number of marker loci in the two files should be the same." )
   return(TRUE)
 }



  return(FALSE)

}









calculateMMt <- function(geno=NULL, availmemGb, ncpu, selected_loci=NA, dim_of_ascii_M=NULL, quiet = TRUE, message=message)
{
 ## internal function to AM
 ## R interface to Rcpp code to calculate M %*% t(M)
 ## Args
 ##      geno        absolute path + file name of binary packed M file
 ##      availmemGb    amount of memory in Gbytes available for creation of MMt
 ##      ncpu    number of cores for matrix operations
 ##      selectedloci an integer vector that gives the column number (0- L-1 ) of the loci that
 ##                   have been selected to act as fixed QTL effects in the model. 
 ##      dim_of_ascii_M    numeric vector with the row, column numbers of M. 
  #------------------------------------------
  # ascii file about to be overwritten
  #------------------------------------------


  if(!file.exists(geno)){
    message(" Error: The binary packed file ", geno, " cannot be found.\n")
    message(" calculateMMt has terminated with errors.")
    return(NULL)
   }
  if(!any(is.na(selected_loci))) selected_loci <- selected_loci-1
  MMt <- calculateMMt_rcpp( f_name_ascii=geno, selected_loci = selected_loci,
                               max_memory_in_Gbytes=availmemGb, num_cores=ncpu, 
                               dims= dim_of_ascii_M, quiet = quiet, message=message) 
  return(MMt)

}  ## end function









calculateMMt_sqrt_and_sqrtinv <- function(MMt=NULL, checkres=TRUE, 
                                           quiet = TRUE , ngpu=0, message=message)
{
  ## internal function to AM
  ## R function for calculating the square root of M * M^t
  ## and the inverse of the square root of MMt
  ## where M * M^t has already been created. 
  ## Using SVD for calculation
  ## Args
  ##  MMt   a matrix of M * M^t
  ##  checkres  when true, the accuracy of the inversion is checked. 

  ## testing that MMt is positive definite
  if(!is.positive.definite(MMt)){
    message(" Error: the matrix multiplication M %*% t(M) is not positive definite. \n")
    message("        This can occur if there are individuals with identical marker \n")
    message("        information. Please remove individuals with identical marker \n")
    message("        information, remembering also to remove their associated phenotype \n")
    message("        information as well. \n")
    message(" Internal function: calculateMMt_sqrt_and_sqrtinv has terminated with errors")
    return(NULL)
  } 
   res <- list()

      MMt.eigen <- eigen(MMt, symmetric=TRUE )
      sqrt_evals <- diag(sqrt(MMt.eigen$values))
      res[["sqrt"]] <- MMt.eigen$vectors %*% sqrt_evals %*% t(MMt.eigen$vectors)
      rm(MMt.eigen, sqrt_evals)
      gc()
      res[["invsqrt"]] <- chol2inv(chol(res[["sqrt"]]))



   if(checkres){
       a <- (res[["sqrt"]] %*% res[["invsqrt"]] )
       if(trunc(sum(diag(a))) != nrow(MMt))
       {
         message(" \n\n\nWARNING: these results may be unstable.\n")
         message(" The sum of the diagonal elements of the square root of M %*% t(M) and its inverse is ", sum(diag(a)), " where \n")
         message("  it should have been ", nrow(MMt), "\n")
         message("  This can occur if the genotype file contains near identical rows and/or columns.  Please check.\n\n")
      

       } 
   }   ## end if(checkres)
   res <- list(sqrt_MMt=res[["sqrt"]], inverse_sqrt_MMt=res[["invsqrt"]] )
  


} ## end function



calculateH <- function(MMt=NULL, varE=NULL, varG=NULL, message=message )
{
  ## internal function to AM
  ## R function for calculating the H variance matrix 
  ## which is
  ##  H = \sigma^2_E I  + \sigma^2_G  MMt
  ## Args:
  ##     MMt  - matrix object for M %*% M^T
  ##     varE  -  numeric value for the residual variance
  ##     varG  -  numeric value for the polygenic variance (\sigma^2_g)
  ##
  ## H matrix is returned. 

  if(!is.numeric(varE)){
    message(" The varE (residual variance) must be numeric.")
    return(NULL)
    }

  if(varE < 0){
    message(" VarE cannot be negative.")
    return(NULL)
    }
  if(!is.numeric(varG)){
    message(" The varG (genotypic variance) must be numeric.")
    return(NULL)
    }
  if(varG < 0){
    message(" VarG cannot be negative.")
    return(NULL)
    }

  if(is.null(MMt)){
    message("MMt cannot be null.")
    return(NULL)
    }
  return( varE * diag(nrow(MMt)) + varG * MMt)


}


calculateP  <- function(H=NULL, X=NULL, ngpu=0, message=message)
{
  ## internal function to AM
  ## R function to calculate P matrix
  ## Args:
  ##       H is the variance matrix
  ##       X is the design matrix supplied by the user
  ##       ngpu for the number of gpu
  ## Returns:
  ##   matrix object P

  if(is.null(H)){
    message(" H must be specified.")
    return(NULL)
  }
  if(is.null(X)){
    message(" A design matrix has not be specified. ")
    return(NULL)
  }

   if(nrow(H) != nrow(X)){
      message(" The number of rows in H and X are not the same.")
    return(NULL)
  }

 Hinv <- chol2inv(chol(H))
 P <- Hinv - Hinv %*% X %*% solve( t(X) %*% Hinv %*% X )  %*% t(X) %*% Hinv

  return(P)

}


calculate_reduced_a <- function(varG=NULL, P=NULL, MMtsqrt=NULL, y=NULL, quiet=TRUE, message=message)
{

  ## internal function to AM
  if( !(nrow(P) ==  length(y))){
    message(" Error:  there is a problem with the  dimensions of  P, and/or the vector y.")
    message("         They should  be of the dimension (n x n), and a vector of length n.")
    message(" The dimensions are: \n")
    message(" dim(P)      = ", dim(P), "\n")
    message(" length(y)   = ", length(y), "\n")
    return(NULL)

  }

 if(is.null(varG)){
   message(" VarG must be specified.")
   return(NULL)
   }

  if(is.null(P)){
   message(" P must be specified")
   return(NULL)
   }


  if(is.null(y)){
   message(" y must be specified")
   return(NULL)
   }

    a <- varG * MMtsqrt %*% P %*% y

return(a)

}






calculate_a_and_vara <- function(geno=NULL, maxmemGb=8, 
                         selectedloci = NA,
                         invMMtsqrt=NULL, transformed_a=NULL, transformed_vara=NULL,
                         quiet = TRUE, message=message)
{
 ## internal function to AM
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



  fnameMt <- geno[["asciifileMt"]]
  dimsMt <- c(geno[["dim_of_ascii_M"]][2], geno[["dim_of_ascii_M"]][1])

  if(!any(is.na(selectedloci))) selectedloci <- selectedloci-1
  calculate_a_and_vara_rcpp(f_name_ascii=fnameMt,
                    selected_loci = selectedloci,
                    inv_MMt_sqrt=invMMtsqrt,  
                    dim_reduced_vara = transformed_vara,
                    max_memory_in_Gbytes=maxmemGb, 
                    dims=dimsMt, 
                    a = transformed_a, 
                    quiet = quiet, message=message)
                    

}


calculate_reduced_vara <- function(X=NULL, varE=NULL, varG=NULL, invMMt=NULL, MMtsqrt=NULL, quiet=TRUE, message=message)
{
 ## internal function to AM
 ## Using var(\hat(a)) = simgaG - Cjj  where Cjj is the component from C^-1 (Henderson's 
 ##   mixed model equations coefficient matrix.   See Verbyla et al. TAG 2007.

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


  vars <- varG * diag(nrow(D1))  - ( D1 + D1 %*% C %*% solve(A - B %*% D1 %*% C) %*% B %*% D1 )

    return(vars )

}



check.inputs <- function(ncpu=NULL, availmemGb=NULL, 
                         file_genotype=NULL, 
                         file_phenotype=NULL )
{
 ## internal function to AM

if(!is.null(ncpu)){
 if(!is.numeric(ncpu)){
   message("Error:  ncpu is not a numeric value. It is of class ", class(ncpu), "It should be the number of cpu.\n")
   return(TRUE)
 }
 if(ncpu < 1){
    message("Error: ncpu cannot be a zero or a negative number. It should be the number of cpu. \n")
    return(TRUE)
 }
    
}

if(!is.null(availmemGb))
{
 if(!is.numeric(availmemGb)){
   message("Error: availmemGb is not a numeric value. It is of class ", class(availmemGb), "It should be the number of gigabytes of RAM available. \n")
   return(TRUE)
 }
 if(availmemGb <= 0){
    message("Error: availmemGb cannot be zero or a negative number.  It should be the number of gigabytes of RAM available. \n")
    return(TRUE)
  } 
}


if(!is.null(file_genotype))
{
  genofile <- fullpath(file_genotype)

  if(!file.exists(genofile)){
    message("Error: Cannot find marker file ", genofile, "\n")
    message("       This could be a problem with the name of the file and/or the location of the file. \n")
    message("       Perhaps specify the full name of the file (i.e. absolute directory path and file name) \n")
    message("       Type help(ReadMarker) and go to the examples section for an example of this. \n") 
    return(TRUE)
  }
}


if(!is.null(file_phenotype))
{ 
  phenofile <- fullpath(file_phenotype)

  if(!file.exists(phenofile)){
    message("Error: Cannot find phenotype file ", phenofile, "\n")
    message("       This could be a problem with the name of the file and/or the location of the file. \n")
    message("       Perhaps specify the full name of the file (i.e. absolute directory path and file name) \n")
    message("       Type help(ReadPheno) and go to the examples section for an example of this. \n") 
    return(TRUE)
  }
}

  return(FALSE)


}




create.ascii  <- function(file_genotype=NULL,  type="text", AA=NULL, AB=NULL, BB=NULL, 
                         availmemGb=8, dim_of_ascii_M=NULL, quiet=TRUE, missing=NULL){
 ## an Rcpp function to create the no-space file of the genotype data M and Mt
 ## from marker data. The marker data may be from an ASCII file or PLINK ped file.
 ## Args
 ## file_genotype    absolute path and file name of genotype file
 ## AA, AB, BB       numeric codes for associated genotypes in marker genotype file
 ## availmemGb     available memory for conversion to reformatted file
 ## dim_of_ascii_M             row, column dimensions of M.  
 ## type            where file type is text or PLINK

 if(.Platform$OS.type == "unix") {
       ##asciiMfile <- paste(dirname(file_genotype), "/", "M.ascii", sep="")
       ##asciiMtfile <- paste(dirname(file_genotype), "/", "Mt.ascii", sep="")
       asciiMfile <- paste(tempdir() , "/", "M.ascii", sep="")
       asciiMtfile <- paste(tempdir()  , "/", "Mt.ascii", sep="")
 } else {
       ##asciiMfile <- paste(dirname(file_genotype), "\\", "M.ascii", sep="")
       ##asciiMtfile <- paste(dirname(file_genotype), "\\", "Mt.ascii", sep="")
       asciiMfile <- paste(tempdir() , "\\", "M.ascii", sep="")
       asciiMtfile <- paste(tempdir() , "\\", "Mt.ascii", sep="")
 }



if (type=="text"){
    ## text genotype file
    if(!is.null(missing)) {
        missing <- as.character(missing)
    } else {
      missing <- "NA"
    }
    it_worked <- createM_ASCII_rcpp(f_name = file_genotype, type=type ,  f_name_ascii = asciiMfile, AA = AA, AB = AB, BB = BB,
               max_memory_in_Gbytes=availmemGb,  dims = dim_of_ascii_M , 
               quiet = quiet, message=message, missing=missing)
    if(!it_worked) #  creation of ASCII file has failed 
       return(FALSE)
 
 
    createMt_ASCII_rcpp(f_name = asciiMfile, f_name_ascii = asciiMtfile,   type=type, 
                  max_memory_in_Gbytes=availmemGb,  dims = dim_of_ascii_M, quiet = quiet, message=message )
} else {
    ## PLINK ped file
    ## using -9 to indicate missing/null genotypes
    ncol  <- dim_of_ascii_M[2] 
    dim_of_ascii_M[2] <- 2*dim_of_ascii_M[2] + 6  ## number of cols in a PLINK file
    it_worked <- createM_ASCII_rcpp(f_name = file_genotype, type=type,  f_name_ascii = asciiMfile, AA ="-9", AB = "-9", BB = "-9",
               max_memory_in_Gbytes=availmemGb,  dims = dim_of_ascii_M , quiet = quiet,   
               message=message, missing="NA")
     if(!it_worked) #  creation of ASCII file has failed 
       return(FALSE)


    dim_of_ascii_M[2] <- ncol ## setting back to number of cols in no-space ASCII file
    createMt_ASCII_rcpp(f_name = asciiMfile, f_name_ascii = asciiMtfile,    type=type, 
                  max_memory_in_Gbytes=availmemGb,  dims = dim_of_ascii_M, quiet = quiet, message=message )

}  ## end if else type

 return(TRUE)

}




extract_geno <- function(fnameM=NULL, colnum=NULL, availmemGb=8, 
                          dim_of_ascii_M=NULL, 
                          selected_locus=NA)
  { ## internal function for AM
    ## Rcpp function to extra a column of genotypes from ascii file M

    selected_locus <- colnum - 1  ## to be consistent with C++'s indexing starting from 0

    genodata <- extract_geno_rcpp(f_name_ascii=fnameM,
                               max_memory_in_Gbytes = availmemGb, 
                              selected_locus=selected_locus, dims=dim_of_ascii_M)

    return(genodata)

  }



constructX <- function(fnameM=NULL, currentX=NULL, loci_indx=NULL, 
                       availmemGb=8, dim_of_ascii_M=NULL,
                        map=NULL)
  {
    ## internal function for AM
    ## R function to construct the design matrix X
    ## Args
    ##   currentX    current model matrix
    ##   loci        the marker loci to be included as fixed QTL effects (additive model)
   
   if(is.na(loci_indx))
   {
     return(currentX)
   } else {
       genodat <- extract_geno(fnameM=fnameM, colnum=loci_indx, 
                           availmemGb=availmemGb, dim_of_ascii_M=dim_of_ascii_M)
      newX <- cbind(currentX, genodat)
      colnames(newX) <- c(colnames(currentX), as.character(map[[1]][loci_indx])) ## adding col names to new X  
      return(newX)
   }
  }





