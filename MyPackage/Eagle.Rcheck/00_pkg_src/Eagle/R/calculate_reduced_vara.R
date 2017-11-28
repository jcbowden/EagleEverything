
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



