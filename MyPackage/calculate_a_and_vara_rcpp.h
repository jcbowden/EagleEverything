//calculate_a_and_vara_rcpp.h

#ifndef calculate_a_and_vara_rcpp_h_INCLUDED   // if x.h hasn't been included yet...
#define calculate_a_and_vara_rcpp_h_INCLUDED   //   #define this so the compiler knows it has been included


Eigen::MatrixXd  ReadBlock(std::string asciifname,
                           long start_row,
                           long numcols,
                           long numrows_in_block);



#endif 
