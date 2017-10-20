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



