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


