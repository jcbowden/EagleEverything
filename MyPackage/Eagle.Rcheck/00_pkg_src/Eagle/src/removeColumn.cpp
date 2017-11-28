// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


const size_t bits_in_int = std::numeric_limits<int>::digits;







void removeColumn(Eigen::MatrixXd& matrix, unsigned long colToRemove)
{
    unsigned long numRows = matrix.rows();
    unsigned long numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}



