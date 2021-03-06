// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


const size_t bits_in_int = std::numeric_limits<int>::digits;





// recode ascii as ascii but with no spaces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool  CreateASCIInospace(std::string fname, std::string asciifname, std::vector<long> dims,
                         std::string  AA,
                         std::string AB,
                         std::string BB,
                         bool  quiet,
                         Rcpp::Function message,
                         std::string missing)
{


std::string
   tmp,
   token,
   line;

 std::ostringstream
      os;

// open marker text  file
std::ifstream fileIN(fname.c_str());

if(!fileIN.good()) {
  message("ERROR: Text file could not be opened with filename  " , fname , "\n" );
  return false;
}
// open ascii file that is to hold  genotype data
 std::ofstream fileOUT(asciifname.c_str(), std::ios::out );
 if (!quiet ){
 message("");
 message(" Reading text File  ");
 message("");
 message(" Loading file ");
 }
long
   number_of_columns,
   counter = 0;



 // initializing input line 

std::vector<char> rowinfile( dims[1] );


 for(long i=0; i < dims[1]; i++)
    rowinfile[i] = '0';



while(getline(fileIN, line ))
{


   Rcpp::Rcout << "\r" << 100.0*counter/dims[0] << "% read of text file.       " << std::flush;

 // Here, BB is coded into 2 
 //       AB is coded into 1, 
 //       AA is coded into 0. 
  std::istringstream streamA(line);
  long i=0;
  number_of_columns = 0;
  while(streamA >> token)
  {
        number_of_columns++;


        if(token == BB){
             rowinfile[i] = '2';
        } else if (token == AB) {
             rowinfile[i] = '1';
        } else if (token == AA) {
             rowinfile[i] = '0';
        } else if (token == missing){
            // setting any missing genotypes to hets
             rowinfile[i] = '1';
        } else {
          if (AB=="NA"){
              message( "\n Marker file contains marker genotypes that are different to AA=" , AA , " BB=" , BB);
              message(" For example , " , token, " in row ", counter+1 );
              message("\n ReadMarker has terminated with errors\n");
              return false;
          } else {
              message( "\n Marker file contains marker genotypes that are different to AA=" , AA , " AB=" , AB , " BB=" , BB);
              message( " For example , " , token , " in row ", counter+1);
              message( "\n ReadMarker has terminated with errors\n");
              return false;
         }
       }  //end if else 
       i++;
  } // end whle streamA


  if (number_of_columns != dims[1] ){
      message("\n");
      message("Error:  Marker text file contains an unequal number of columns per row.  ");
      message("        The error has occurred at row " , counter+1 , " which contains " , number_of_columns , " but ");
      message("        it should contain " , dims[1] , " columns of data. ");
      message("\n");
      message(" ReadMarkerData has terminated with errors");
      return false;
  }
  for(long ii=0; ii< number_of_columns; ii++){
     fileOUT << rowinfile[ii];
  }
  fileOUT << "\n";
  counter++;

 }  // end while getline



  // write out a few lines of the file 
     std::ifstream fileIN_backtobeginning(fname.c_str());
     counter = 0;

    int nrowsp =  5;
    int ncolsp = 12;
    if(dims[0] < 5)
    nrowsp = dims[0];
    if(dims[1] < 12)
    ncolsp = dims[1];

    message(" First ", nrowsp, " lines and ", ncolsp, " columns of the marker text  file. ");

    std::string rowline;
    while(getline(fileIN_backtobeginning, line ) && counter < nrowsp)
    {
       std::ostringstream oss;
       std::istringstream streamA(line);
       for(int i=0; i < ncolsp ; i++){
           streamA >> tmp;
           oss << tmp << " " ;
        }
        std::string rowline = oss.str();
        message(rowline);
        counter++;
      }  // end  while(getline(fileIN, line ))

// close files
fileIN.close();
fileOUT.close();

  return true;
}



