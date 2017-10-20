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


