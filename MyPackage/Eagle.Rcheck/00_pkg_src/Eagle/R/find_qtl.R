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



