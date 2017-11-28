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


