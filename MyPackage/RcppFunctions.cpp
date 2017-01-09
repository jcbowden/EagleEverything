// This software is distributed under the GNU General Public License.
//
//This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 


// Author:   Andrew W. George
// Purpose: to calculate M %*% t(M) when M may not fit into memory
// Outline: 
//          1. read data from PLINK or text file
//          2. convert genotypes into their binary values.
//          3. pack binary values into unsigned long int (could be 32 bits or 64 bits 
//             depending upon the system.
//          4. write packed longs to a new file in binary format.
//          5. read blocks of binarys to form submatrices of M.
//          6. Perform M %*% t(M) as a block multiplication. 
//
#define EIGEN_USE_MKL_ALL

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <R.h>

#include <omp.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <bitset>
#include <string>

#include <ctime>

//    #include <magma.h>



using namespace std;
using namespace Rcpp;

using Eigen::MatrixXi;
using Eigen::MatrixXd;  
using Eigen::Lower;
using Eigen::Map;   // maps rather than copies

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


const size_t bits_in_double = std::numeric_limits<long double>::digits;
const size_t bits_in_ulong = std::numeric_limits<unsigned long int>::digits;
const size_t bits_in_int = std::numeric_limits<int>::digits;






// check genotypes in file are correct numeric values AA, AB, BB
// [[Rcpp::export]]
void  checkGenotypes(CharacterVector f_name,
                     int AA,
                     int AB,
                     int BB,
                     bool csv)
{
std::string 
     fname = Rcpp::as<std::string>(f_name),
     token,
     line;
 int 
    genoval=-1,
    linenum=0;

 ostringstream 
      os;

 char 
   sep = ' ';
 if(csv) sep = ',';


 // open file and check for its existence. 
 std::ifstream fileIN(fname.c_str());
 if(!fileIN.good()) {
      os << "\n\nERROR: Could not open  " << fname << "\n\n" << std::endl;
      Rcpp::stop(os.str() );
 }

 // Determine number of rows in file
 Rprintf("\n\n Checking genotype file for incorrect genotypes ... ");
 while(fileIN.good()){
      while(getline(fileIN, line)){
           Rprintf(".");
           linenum++;
           istringstream streamA(line);
          //  while(streamA >> genoval)
           // while(getline(streamA,  genoval, sep))
            while(getline(streamA,  token, sep)){
                genoval = atoi(token.c_str());
           //  if(genoval !=0 & genoval !=1 & genoval != 2){
             if(genoval !=AA & genoval !=AB & genoval != BB){
               os << "\n\nERROR: File " << fname << " contains genotypes other than " << AA << "," << 
                        AB << ", and " << BB << " For example genotype " << genoval << " has been found on line.  " << linenum << "\n\n";
               Rcpp::stop(os.str() );
             }
           }
      }
 }


}





//get number of rows and columns in marker file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// [[Rcpp::export]]
std::vector<long>   getRowColumn(std::string fname, 
                                 bool csv)
{
  // Purpose:  to open the marker file where the marker data are kept.
  //           An error will be produced if the file cannot be found.
  //           I am assuming no row or column names
 int 
   genoval;

 std::string
   line;

 ostringstream 
      os;


 std::vector<long> dimen(2,0)  ;  // dim[0] row number
                               // dim[1] col number 

 char 
    sep = ' ';
 if(csv)
     sep = ',';

 // open file and check for its existence. 
 std::ifstream fileIN(fname.c_str());
 if(!fileIN.good()) {
      os << "\n\n ERROR: Could not open  " << fname << "\n\n" << std::endl;
      Rcpp::stop(os.str() );
 }


 // Determine number of rows in file
 while(fileIN.good()){
      while(getline(fileIN, line )){
         dimen[0]++;
      }
 }


 // Determine number of columns in file
 fileIN.clear(); // returns to beginning of line
 fileIN.seekg(0, ios::beg);

 getline(fileIN, line );
 istringstream streamA(line);

  string 
    token ;

 while(streamA >> token)
   dimen[1]++;
 // while(getline(streamA , token , sep)){
 //          dimen[1]++;
//  }
 fileIN.close();
 if(fileIN.bad())
 {
    os << "\n\nERROR:  There was a problem with reading the marker file - possibly strange ASCII characters.\n\n";
    Rcpp::stop(os.str() );
}
return dimen;

}






// recode PLINK as ASCII with no spaces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void  CreateASCIInospace_PLINK(std::string fname, std::string asciifname, std::vector<long> dims,
                         bool quiet)
{
long 
   colindx = 0; 

int
  n_of_cols_in_geno = (dims[1] -6)/2.0;

int
  genovec[n_of_cols_in_geno];

char
   alleles [ 2 ][ n_of_cols_in_geno ];  // holds alleles  


std::vector<char>
     rowvec( dims[1] - 6 );  // holds allelic information from PLINK file


std::string
   tmp,
   token,
   line;

char 
   sep = ' ';


 ostringstream 
      os;




// open PLINK ped  file
std::ifstream fileIN(fname.c_str());
if(!fileIN.good()) {
  Rcpp::Rcout << "\n\nERROR: PLINK ped file could not be opened with filename  " << fname <<  std::endl;
  os << "\n\nERROR: ReadMarkerData has terminated with errors.  " << fname << "\n\n" << std::endl;
  Rcpp::stop(os.str() );
}

// open ascii file that is to hold no-spaces genotype data
std::ofstream fileOUT(asciifname.c_str(), ios::out );
long  counter = 0;
long  number_of_columns;
Rcout << "\n Reading marker file " ;

// initializing input line 
std::string rowinfile(n_of_cols_in_geno, '0'); // s == "000000"




while(getline(fileIN, line ))
{
    Rcout << " Processing line  " << counter << endl;

    istringstream streamA(line);


    for(long i=0; i < dims[1] ; i++){
          // assign allelic info to rowvec ignoring first 6 columns of input
         if(i <= 5){
                 streamA >> tmp;
          } else {
                 streamA >> rowvec[i-6];
          }
    }  // end  for(long i=0; i < dims[1] ; i++)


    // initialize alleles structure to first row of PLINK info
    if (counter == 0) {
         for(long i=0; i < n_of_cols_in_geno ; i++){
              alleles[ 0 ][ i ] =  rowvec[ (2*i ) ];
              alleles[ 1 ][ i ] =  rowvec[ (2*i + 1) ];
             //   Rcout << "alleles "  << alleles[0][i] << alleles[1][i] << endl;
         }
    }

    // turn allelic info from PLINK into genotype 0,1,2 data
    // also do some checks for more than 2 alleles, and 0 and - for missing data
    for(long i=0; i < n_of_cols_in_geno; i++){
          // Checking for missing allelic information in PLINK file
         if( rowvec[ (2*i ) ] == '0' ||  rowvec[ (2*i + 1) ] == '0' || rowvec[ (2*i ) ] == '-' ||  rowvec[ (2*i + 1) ] == '-'){
                 Rcpp::Rcout << std::endl;
                 Rcpp::Rcout << std::endl;
                 Rcpp::Rcout << "Error:  PLINK file cannot contain missing alleles (i.e. 0 or - ) " << std::endl;
                 Rcpp::Rcout << "        Please impute missing marker information before running AMplus." << std::endl;
                 Rcpp::Rcout << "        The error has occurred at snp locus " << i << " for individual " << counter+1 << std::endl;
                 Rcpp::Rcout << std::endl;
                 Rcpp::Rcout << std::endl;
                 os << " ReadMarkerData has terminated with errors\n" << std::endl;
                 Rcpp::stop(os.str() );
        }   

        // Check if allele has been seen before in allele file. 
        // If so, make sure alleles doesn't already  contain two alleles - otherwise generate error message
        for(int j = 1; j >= 0; --j){ // looping over the two alleles with indexes 0 and 1
               if (rowvec[ (2*i + j) ] != alleles[ 0 ][ i ] && rowvec[ (2*i + j) ] != alleles[ 1 ][ i ]){
                     if (alleles[ 0 ][ i ] == alleles[ 1 ][ i ] ){
                           // this is okay. alleles only contains a single allele at the moment. Re-initialise alleles
                           alleles[ 1 ][ i ] = rowvec[ (2*i + j) ];
                     } else {
                           // Error - we have more than two alleles segregating at a locus
                           Rcpp::Rcout << std::endl;
                           Rcpp::Rcout << std::endl;
                           Rcpp::Rcout << "Error:  PLINK file cannot contain more than two alleles at a locus."  << std::endl;
                           Rcpp::Rcout << "        The error has occurred at snp locus " << i << " for individual " << counter+1 << std::endl;
                           Rcpp::Rcout << std::endl;
                           Rcpp::Rcout << std::endl;
                           os << " ReadMarkerData has terminated with errors\n" << std::endl;
                           Rcpp::stop(os.str() );
                    } // end inner if else
               }  // end if
       }  // end for(int j = 1; j >= 0; --j)


    // set genovec
    if (rowvec[ (2*i + 1) ] !=   rowvec[ (2*i) ] ){
         genovec[i] = 1 ;  // AB
    } else {
         if (rowvec[ (2*i ) ] == alleles[ 0 ][ i ] ){  // matches first allele
                genovec[i] = 0;  // AA
          }  else {
                 genovec[i] = 2;  // BB
         }
    }  // end outer if else rowvec
  //  fileOUT << genovec[i];

  } // end  for(long i=0; i < n_of_cols_in_geno; i++)

 counter++;
 fileOUT << genovec;
 fileOUT << "\n";

}  // end while(getline(fileIN, line ))
Rcout << "\n\n" << endl;
// write out a few lines of the file if quiet
// open PLINK ped  file
std::ifstream fileIN_backtobeginning(fname.c_str());
counter = 0;
Rcout << " First 5 lines and 12 columns of the PLINK ped  file. " << endl;
while(getline(fileIN_backtobeginning, line ) && counter < 5)
{
       Rcpp::Rcout << " " ;
       istringstream streamB(line);
       for(int i=0; i < 12 ; i++){
           streamB >> tmp;
           Rcpp::Rcout << tmp << " " ;
        }
        Rcpp::Rcout << std::endl;
        counter++;
}  // end  while(getline(fileIN, line ))





// close files
fileIN.close();
fileOUT.close();


}











// recode ascii as ascii but with no spaces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void  CreateASCIInospace(std::string fname, std::string asciifname, std::vector<long> dims,
                         std::string  AA, 
                         std::string AB, 
                         std::string BB,
                         bool csv, 
                         bool quiet)
{
long 
   colindx = 0;

int 
   genovec[dims[1]];


 std::string 
     rowvec[dims[1]]; // holds entire row worth of genotypes from ascii file
 

std::string
   tmp,
   token,
   line;
//    rowinfile;

char 
   sep = ' ';
if(csv) 
   sep = ',';


 ostringstream 
      os;


// open marker text  file
std::ifstream fileIN(fname.c_str());

if(!fileIN.good()) {
  os << "\n\nERROR: Text file could not be opened with filename  " << fname << "\n\n" << std::endl;
  Rcpp::stop(os.str() );
}

// open ascii file that is to hold  genotype data
 std::ofstream fileOUT(asciifname.c_str(), ios::out );
 if (quiet > 0){
 Rcpp::Rcout << " " << std::endl;
 Rcpp::Rcout << " Reading text File  " << std::endl;
 Rcpp::Rcout << " " << std::endl;
 Rcpp::Rcout << " Loading file .";
 }
long 
   number_of_columns, 
   counter = 0;


 // initializing input line 
 std::string rowinfile(dims[1], '0'); // s == "000000"


while(getline(fileIN, line ))
{
  Rcout << " Processing line  " << counter << endl;


 // Here, BB is coded into 2 
 //       AB is coded into 1, 
 //       AA is coded into 0. 
  istringstream streamA(line);
  long i=0;
  number_of_columns = 0;
  while(streamA >> token)
  {
 //    if(quiet > 0)
        number_of_columns++;


        if(token == BB){
             rowinfile[i] = '2';
        } else if (token == AB) {
             rowinfile[i] = '1';
        } else if (token == AA) {
             rowinfile[i] = '0';
        } else {
          if (AB=="NA"){
              Rcpp::Rcout << "Error: Marker text file contains marker genotypes that are different to " << AA << " " << BB << endl;
              Rcpp::Rcout << "       For example, " << rowvec[i] << endl;
              os << "ReadMarker has terminated with errors\n\n"; 
              Rcpp::stop(os.str() );
          } else {
              Rcpp::Rcout << "Error: Marker text file contains marker genotypes that are different to " << AA << " " << AB << " " << BB << endl;
              Rcpp::Rcout << "       For example, " << rowvec[i] << endl;
              os << "ReadMarker has terminated with errors\n\n"; 
              Rcpp::stop(os.str() );
         }
       }  //end if else 
       i++;
  } // end whle streamA

  if (quiet>0){
        Rcout << " Number of columns in line " << counter+1 << " is " << number_of_columns << std::endl;

        if (number_of_columns != dims[1] ){
             Rcpp::Rcout << std::endl;
             Rcpp::Rcout << std::endl;
             Rcpp::Rcout << "Error:  Marker text file contains an unequal number of columns per row.  " << std::endl;
             Rcpp::Rcout << "        The error has occurred at row " << counter+1 << " which contains " << number_of_columns << " but " << endl;
             Rcpp::Rcout << "        it should contain " << dims[1] << " columns of data. "  << std::endl;
             Rcpp::Rcout << std::endl;
             Rcpp::Rcout << std::endl;
             os << " ReadMarkerData has terminated with errors\n" << std::endl;
             Rcpp::stop(os.str() );
       }  // end if number_of_columns
  } // end if quiet

  fileOUT << rowinfile;
  fileOUT << "\n";
  counter++;

 }  // end while getline

 Rcout << " FINISHED WRITING M.ascii FILE !!!!!"   << endl; 

 if (quiet > 0) 
     Rcpp::Rcout << "\n" << std::endl;


  // write out a few lines of the file if quiet
  if(quiet > 0){
     // open PLINK ped  file
     std::ifstream fileIN(fname.c_str());
     counter = 0;
     Rcout << " First 5 lines and 12 columns of the text file. " << endl;
     while(getline(fileIN, line ) && counter < 5)
     {
       istringstream streamA(line);
       for(int i=0; i < 12 ; i++){
           streamA >> tmp;
           Rcpp::Rcout << tmp << " " ;
        }
        Rcpp::Rcout << std::endl;
        counter++;
      }  // end  while(getline(fileIN, line ))
  } // end if(quiet)


// close files
fileIN.close();
 fileOUT.close();
// fclose(pfileOUT);

}










Eigen::MatrixXd  ReadBlock(std::string asciifname, 
                           long start_row,
                           long numcols,
                           long numrows_in_block)

{
 // reads in data from ASCII file 
 // to form M Eign double matrix 
ostringstream
      os;
std::string
   line;



long 
  coli = 0, 
  rowi = 0;

long 
    igeno;

double 
     geno;

Eigen::MatrixXd
      M(numrows_in_block, numcols) ;


// Open no-space ASCII file
   std::ifstream fileIN(asciifname.c_str(), ios::in );

    if(!fileIN.good()) {
      os << "ERROR: Could not open  " << asciifname << std::endl;
      Rcpp::stop(os.str() );
     }

   for(long rr=0; rr < (start_row + numrows_in_block) ; rr++){
      // read a line of data from ASCII file
      getline(fileIN, line);
      if(rr >= start_row){
          istringstream streamA(line);
          for(long ii=0; ii < numcols  ; ii++){
            M(rowi, ii)  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
          }
          rowi++;
      } // end if rr
   } // end for(rr



// Close the binary file
   fileIN.close();


 return M;

}








// [[Rcpp::export]]
void  createMt_ASCII_rcpp(CharacterVector f_name, CharacterVector f_name_ascii, 
                              double  max_memory_in_Gbytes,  std::vector <long> dims,
                              bool quiet )
{

// read data from M.ascii that has already been created and transpose this file

std::string
   token, 
   line;

ostringstream
      os;

long
 n_of_cols_to_be_read;

int
   genoval;

std::string
     fname = Rcpp::as<std::string>(f_name),
     fnameascii = Rcpp::as<std::string>(f_name_ascii);

std::string 
  rowvec[dims[1]];

char 
   sep = ' ';


double 
   max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;


// Calculate number of columns that can be read in as a block with XGb of
// memory. 
   // how much memory will be needed to  store M, takes it transpose M.transpose, and 
   // store its answer in Mt (+ .5 for a buffer)
double mem_bytes = 3.5 * dims[0] * dims[1] * (bits_in_int/8);  // assumes a 64 bit system


 // open  files
std::ofstream fileOUT(fnameascii.c_str(), ios::out);
std::ifstream fileIN(fname.c_str());

//-----------------------------------------------------------------------
//  Two situations
//   1.  memory X is sufficient to read all data into memory and transpose
//   2.  memory X is insufficient to read all data into memory. 
//------------------------------------------------------------------------



if(mem_bytes < max_mem_in_bytes){
     // Situation 1
     //-------------

     // open ASCII file and check for its existence. 
     if(!fileIN.good()) {
        os << "\n\nERROR: Could not open  " << fname << "\n\n" << std::endl;
        Rcpp::stop(os.str() );

     }

    // create matrix structure to hold genotype data
    Eigen::MatrixXi 
          M(dims[0], dims[1]) ;


    // reset position in data file
    fileIN.clear();
    fileIN.seekg(0, ios::beg);

   // read values into matrix
   long rowi=0;
   Rcout << " Reading in M.ascii into M matrix " << endl;
   while(getline(fileIN, line ))
   {
       // read a line of data from ASCII file
       for(long coli=0; coli < dims[1]; coli++){
           M(rowi, coli)  = line[coli] - '0'; // trick to removes ASCII character offset for numbers
       }
       rowi++;
   }  // end while getline

  // take transose of matrix M
  Rcout << " Taking transose " << endl;
  MatrixXi Mt = M.transpose();
  
  // write out contents fo Mt to file (no spaces)a
  Rcout << " WRiting out contents of Mt to file " << endl;
  std::string rowinfile(Mt.cols(), '0');  // initialising string with 0's
  for(long rowi=0; rowi<Mt.rows(); rowi++){
     for(long coli=0; coli<Mt.cols(); coli++){
         // fileOUT << Mt(rowi, coli);
        rowinfile[coli] =  Mt(rowi, coli) + '0'; // forming string row before writing to file
     } 
     fileOUT << rowinfile; // writing entire row of data
     fileOUT << "\n";
  }
 Rcout << " Finished ... " << endl;
 
 fileIN.close();
 fileOUT.close();

} else {
   //  Situation 2 
   //  Block approach needed due to lack of memory

   if (quiet > 0){
        Rcpp::Rcout << " A block transpose is being performed due to lack of memory.  "  << std::endl;
        Rcpp::Rcout << " Memory parameter workingmemGb is set to " << max_memory_in_Gbytes << "Gbytes" << std::endl;
        Rcpp::Rcout << " If possible, increase workingmemGb parameter. " << std::endl;
    }

    // Calculate number of columns that can be read into available memory
    n_of_cols_to_be_read = max_mem_in_bytes * 1.0 / (3.5 * dims[0] * (bits_in_int/8.0)); //64 bit system
    Rcout << n_of_cols_to_be_read << endl;
    // Calculate number of blocks needed
    long n_blocks = dims[1]/n_of_cols_to_be_read;
    if (dims[1] % n_of_cols_to_be_read != 0)
         n_blocks++;
     Rcout << " number of blocks is " << n_blocks << endl;
    if (quiet > 0)  
         Rcpp::Rcout  << " Block Tranpose of ASCII genotype file beginning ... " << std::endl;

    // Block read and transpose - requires n_blocks passes through the 
    // ASCII input file which could be slow if file is large and memory low
    for(long b=0; b < n_blocks; b++){
         if (quiet > 0) 
              Rcpp::Rcout << " Processing block ... " << b << " of a total number of blocks of " << n_blocks << std::endl;


         long
              start_val = b * n_of_cols_to_be_read,
              end_val   = (b+1) * n_of_cols_to_be_read;

         if (end_val > dims[1])
              end_val = dims[1];

         long ncols = end_val - start_val ;
         MatrixXi
              M(dims[0], ncols );


         // open ASCII file and check for its existence. 
         if(!fileIN.good()) {
              os << "ERROR: Could not open  " << fname << std::endl;
              Rcpp::stop(os.str() );
         }
         long counter = 0;
         if (quiet > 0) {
              Rcpp::Rcout << std::endl;
              Rcpp::Rcout << std::endl;
         }
         for(long rowi=0; rowi<dims[0]; rowi++){ 
               // read a line of data from ASCII file
              getline(fileIN, line);
              istringstream streamA(line);
  
             long coli=0 ;
             for(long ii=start_val; ii < end_val ; ii++){
                M(rowi, coli)  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
                coli++;
             }
        } // end for(rowi=0; rowi<dims[0]; rowi++)
       // tranpose M


       MatrixXi Mt = M.transpose();


      // write out contents fo Mt to file (no spaces)a
      std::string rowinfile(Mt.cols(), '0');  // initialising string with 0's
      for(long rowi=0; rowi<Mt.rows(); rowi++){
         for(long coli=0; coli<Mt.cols(); coli++){
             // fileOUT << Mt(rowi, coli);
            rowinfile[coli] =  Mt(rowi, coli) + '0'; // forming string row before writing to file
         }
         fileOUT << rowinfile; // writing entire row of data
         fileOUT << "\n";
      }




      // return to the beginning of the input file
      fileIN.clear();
      fileIN.seekg(0, ios::beg);

     }     // end for blocks

fileIN.close();
fileOUT.close();


}  // end if else situation 





}  // end function 







//--------------------------------------------
// Calculation of transformed blup a values
//--------------------------------------------
// [[Rcpp::export]]
MatrixXd calculate_reduced_a_rcpp ( CharacterVector f_name_ascii, double varG, 
                                           Map<MatrixXd> P,
                                           Map<MatrixXd>  y,
                                           double max_memory_in_Gbytes,  
                                           std::vector <long> dims,
                                           Rcpp::NumericVector  selected_loci,
                                           bool quiet)
{
  // function to calculate the BLUPs for the dimension reduced model. 
  // It is being performed in Rcpp because it makes use of Mt. 
  // Args
  // f_name_ascii    path + file name of Mt.bin
  // varG          variance of polygenic component
  // P             calculate in R
  // y             response/trait  but read in as a row matrix
  // max_memory_in_Gbytes  working memory in Gbytes
  // dims          dimension (row, column), of M.

std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

Eigen::MatrixXd
      ar(dims[1],1);  // column vector

ostringstream
      os;



const size_t bits_in_double = std::numeric_limits<double>::digits;


   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt))
double mem_bytes_needed =   ( dims[0]*dims[1] + dims[0]*dims[0] + dims[0] ) *  ( sizeof(double)/( 1000000000));

if (!quiet){
    // Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
    Rprintf("Max memory (Gbytes) available is: %f \n", max_memory_in_Gbytes);
}

if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory

   Eigen::MatrixXd
                   Mt;


  if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(long ii=0; ii < selected_loci.size() ; ii++)
          Mt.row(selected_loci(ii) ).setZero();
   }

   Rcout << " calculate_reduced_a_rcpp " << " about to read in Mt " << endl;
   Mt = ReadBlock(fnamebin, 0, dims[0], dims[1]);
   Rcout << " read ... " << endl;
   // ar  =    varG * Mt *  P   * y ;
 //    std::clock_t    start;
 //    start = std::clock();
   Rcout << " ar  =    varG * Mt *  P   * y " << endl;
    ar  =     P   * y ;
    ar  =    Mt * ar;
    ar  =    varG * ar;
//    Rcout << "Time1: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
} else {

      // calculation being processed in block form
      Rprintf(" Note:  Increasing workingmemGb would improve performance... \n");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algrebra but it is comes to the following
      long num_rows_in_block = (max_memory_in_Gbytes * 1000000000.0/sizeof(double) - dims[0] * dims[0] - dims[0])/dims[0] ;

    if (num_rows_in_block < 0){
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << "Error:  workingmemGb is set to " << max_memory_in_Gbytes << std::endl;
        Rcpp::Rcout << "        Cannot even read in a single row of data into memory." << std::endl;
        Rcpp::Rcout << "        Please increase workingmemGb for this data set." << std::endl;
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << std::endl;
        os << " multiple_locus_am has terminated with errors\n" << std::endl;
         Rcpp::stop(os.str() );

      }


      // blockwise multiplication

      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;

      if (quiet){
      Rprintf(" Maximum memory has been set to %f Gb\n", max_memory_in_Gbytes);
      Rprintf(" Block multiplication necessary. \n");
      Rprintf(" Number of blocks needing block multiplication is ... % d \n", num_blocks);
      } 
      for(long i=0; i < num_blocks; i++){
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[1])
            num_rows_in_block1 = dims[1] - start_row1;

          Eigen::MatrixXd
                  Mt;
          Mt = ReadBlock(fnamebin, start_row1, dims[0], num_rows_in_block1) ;

         Eigen::MatrixXd
             ar_tmp;

         if(!R_IsNA(selected_loci(0))){
         // setting columns (or row when Mt) to 0
            for(long ii=0; ii < selected_loci.size() ; ii++)
            {
            // since we are now dealing with Mt, and blocking on columns, 
            // because columns are rows in Mt, then we have to be careful
            // that we do not select loci outside the block bounds. Also 
            // the values have to be adjusted based on the block number
                if(selected_loci(ii) >= start_row1 & selected_loci(ii) < start_row1 + num_rows_in_block1 )
                {   // selected loci index is in block 
                long block_selected_loci = selected_loci(ii) - start_row1;
                Mt.row(block_selected_loci).setZero();
                }
             }   
         }

         // ar_tmp  =  varG * Mt *  P  * y ;
          ar_tmp = P * y;
          ar_tmp = Mt * ar_tmp;
          ar_tmp = varG * ar_tmp;
          


         // assign block vector results to final vector (ar) of results
         long  counter = 0;
         for(long j=start_row1; j < start_row1 + num_rows_in_block1; j++){
              ar(j,0) = ar_tmp(counter,0);
              counter++;
         }

       if (quiet)  Rcpp::Rcout << "block done ... " << std::endl;
      } // end for long




}  // end if mem_bytes_needed

  return(ar);

} // end function 




// internal function to remove a row from a dynamic matrix
void removeRow(MatrixXd& matrix, unsigned long rowToRemove)
{
    unsigned long numRows = matrix.rows()-1;
    unsigned long numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}



void removeColumn(Eigen::MatrixXd& matrix, unsigned long colToRemove)
{
    unsigned long numRows = matrix.rows();
    unsigned long numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}






// ------------------------------------------------------
//    Calculation of untransformed BLUP a values 
// ------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List   calculate_a_and_vara_rcpp(  CharacterVector f_name_ascii,  
                                    Rcpp::NumericVector  selected_loci,
                                    Map<MatrixXd> inv_MMt_sqrt,
                                    Map<MatrixXd> dim_reduced_vara,
                                    double  max_memory_in_Gbytes,  
                                    std::vector <long> dims,
                                    Eigen::VectorXd  a,
                                    bool quiet,
                                    Rcpp::NumericVector indxNA)
{
// Purpose: to calculate the untransformed BLUP (a) values from the 
//          dimension reduced BLUP value estimates. 
//          It is neccessary to have a block multiplication form of this function. 
//          Also, since the matrix multiplications are reliant upon the BLAS library, only 
//          double precision matrix multiplication is possible. This means, the Mt matrix must 
//          be converted into a douple precision matrix which has a large memory cost.  
// Note:
//      1. dims is the row, column dimension of the Mt matrix
//      2. when indxNA not NA, then need to adjust dimenions of Mt by removing cols.


ostringstream
      os;

std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

 Eigen::MatrixXd
       ans(dims[0],1);

Eigen::MatrixXd
             ans_tmp,
             var_ans_tmp(dims[0] , dims[1]);



Eigen::MatrixXd
    var_ans = Eigen::MatrixXd(dims[0],1);

const size_t bits_in_double = std::numeric_limits<double>::digits;
const size_t bits_in_integer = std::numeric_limits<int>::digits;



   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt)%*%M)
//AWG  double mem_bytes_needed =   ( dims[0]   +  2*dims[1]   + 1 ) *  (dims[1] * sizeof(double) /( 1000000000));
 double mem_bytes_needed =   ( 4   *dims[1]  *  dims[0] * sizeof(double))/1000000000;

if (!quiet){
//   Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
   Rprintf("Max memory (Gbytes) available is: %f \n", max_memory_in_Gbytes);
}





if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory


    Eigen::MatrixXd Mt = ReadBlock(fnamebin, 0, dims[1], dims[0]);

  // removing columns that correspond to individuals with no 
  // trait data
//   if(!R_IsNA(indxNA(0)))
   if(indxNA.size()!=0){
     for (long ii=0; ii < indxNA.size(); ii++){
        removeColumn(Mt, indxNA(ii) );
     }
  }



   if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(long ii=0; ii < selected_loci.size() ; ii++){
           Mt.row(selected_loci(ii)).setZero();
    }
   }




//   AWGans =    Mtd *  inv_MMt_sqrt  * a ;
//   AWGans =    Mt.cast<double>() *  inv_MMt_sqrt  * a ;
std::clock_t    start;

//   start = std::clock();
    Eigen::MatrixXd  ans_part1 = inv_MMt_sqrt * a;
    ans.noalias() =   Mt  * ans_part1; 

//   Rcout << "Time2 Mtd *  inv_MMt_sqrt  * a : " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;




// AWG    ans_part1.resize(0,0);  // erase matrix
   //  ans =    Mt *  inv_MMt_sqrt  * a ;
 //  Rprintf(" finished untransfomred BLUP values \n");



  // calculate untransformed variances of BLUP values
//  Eigen::MatrixXd var_ans_tmp_part1 =  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;

//     start = std::clock();
    Eigen::MatrixXd var_ans_tmp_part1 =   dim_reduced_vara * inv_MMt_sqrt;
    var_ans_tmp_part1 = inv_MMt_sqrt * var_ans_tmp_part1;
//     Rcout << "Time3 var_ans_tmp_part1 =  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt : " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
//  Eigen::MatrixXd var_ans_tmp_part1 =  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;a

    var_ans_tmp  =  Mt  *  var_ans_tmp_part1;

  var_ans_tmp_part1.resize(0,0);  // erase matrix 

  // Added 26 April
  long i;
  #pragma omp parallel for shared(var_ans, var_ans_tmp, Mt)  private(i) schedule(static)
  for(i=0; i< dims[0]; i++){
           var_ans(i,0) =   var_ans_tmp.row(i)   * (Mt.row(i)).transpose() ;
  }



} else {
    //  -----------------------------------------
    //       BLOCK WISE UPDATE
    //  -----------------------------------------

      // ans.resize(dims[0],1);   //added AWG 12/03/16 in a bid to improve GPU performance

      // calculation being processed in block form
      Rprintf(" Increasing maxmemGb would improve performance... \n");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algrebra but it is comes to the following
      // Strickly, 2 should be used here but I want some extra memory to play with 
      long num_rows_in_block =  max_memory_in_Gbytes * ( 1000000000) /
                             ( 4  *dims[1] *  sizeof(double) ) ;


      if (num_rows_in_block < 0){
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << "Error:  workingmemGb is set to " << max_memory_in_Gbytes << std::endl;
        Rcpp::Rcout << "        Cannot even read in a single row of data into memory." << std::endl;
        Rcpp::Rcout << "        Please increase workingmemGb for this data set." << std::endl;
        Rcpp::Rcout << std::endl;
        Rcpp::Rcout << std::endl;
        os << " multiple_locus_am has terminated with errors\n" << std::endl;
         Rcpp::stop(os.str() );

      }
     


      // blockwise multiplication

      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;
      if (quiet){
      Rprintf(" Maximum memory has been set to %f Gb\n", max_memory_in_Gbytes);
      Rprintf(" Block multiplication necessary. \n");
      Rprintf(" Number of blocks needing block multiplication is ... % d \n", num_blocks);
      } 
      for(long i=0; i < num_blocks; i++){
         Rcpp::Rcout << "Performing block iteration ... " << i << endl;
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[0])
            num_rows_in_block1 = dims[0] - start_row1;


         Eigen::MatrixXd Mt = ReadBlock(fnamebin, start_row1, dims[1], num_rows_in_block1) ;
         // Rcout << Mt.rows() << endl;
         // Rcout << Mt.cols() << endl;
         // removing columns that correspond to individuals with no 
         // trait data
         //  if(!R_IsNA(indxNA(0)))
         if(indxNA.size() != 0 ){
               Rcpp::Rcout << " Removing columns don't know why though ... " << endl;
               for (long ii=0; ii < indxNA.size(); ii++){
                     removeColumn(Mt, indxNA(ii) );
               }
          }



        Eigen::MatrixXd
              vt1,
              ans_tmp1;
          
        Eigen::MatrixXd   var_ans_tmp(num_rows_in_block1,1);

            if(!R_IsNA(selected_loci(0))){
             Rcout << " in if(!R_IsNA(selected_loci(0))) " << endl;
            // setting columns (or row when Mt) to 0
               for(long ii=0; ii < selected_loci.size() ; ii++)
               {
               // since we are now dealing with Mt, and blocking on columns, 
               // because columns are rows in Mt, then we have to be careful
               // that we do not select loci outside the block bounds. Also 
               // the values have to be adjusted based on the block number
                   if(selected_loci(ii) >= start_row1 & selected_loci(ii) < start_row1 + num_rows_in_block1 )
                   {   // selected loci index is in block 
                   long block_selected_loci = selected_loci(ii) - start_row1;
                   Mt.row(block_selected_loci).setZero();
                   }
                }   
            }
           //  ans_tmp  =  Mtd *  inv_MMt_sqrt  * a ;
            Rcout << " 1 ans_tmp.noalias()  =   inv_MMt_sqrt  * a " << endl;
             ans_tmp.noalias()  =   inv_MMt_sqrt  * a ;
         Rcout << "ns_tmp = Mt * ans_tmp " << endl;
             ans_tmp = Mt * ans_tmp;

            // variance calculation
            // vt.noalias() =  Mtd *  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;
Rcout << "vt1.noalias() =  dim_reduced_vara * inv_MMt_sqrt " << endl;
             vt1.noalias() =  dim_reduced_vara * inv_MMt_sqrt;
     Rcout << " vt1           =  inv_MMt_sqrt * vt1 " << endl;
             vt1           =  inv_MMt_sqrt * vt1;



        // performing quadratic form, remembering only diag elements are needed for variances. 
//        Rcout << " Performing ... Mt.row(j) * vt1 * ((Mt.row(j)).transpose())  " << endl;
//        for(long j=0; j < num_rows_in_block1; j++){
//           var_ans_tmp(j,0) =  Mt.row(j) * vt1 * ((Mt.row(j)).transpose()) ;
//        }
//        Rcout << "end of computing variances ... " << endl;
            // vt.noalias() =  Mt *  vt1;
           Eigen::MatrixXd vt;
             Rcout << " vt.noalias()  =  Mt *  vt1 " << endl; 
              vt.noalias()  =  Mt *  vt1;






   //    var_ans_tmp(j,0)  =   vt.row(j)  * ((Mt.row(j)).transpose()) ;
           // Added 26 April
            #pragma omp parallel for
            for(long j=0; j < num_rows_in_block1; j++){
                      var_ans_tmp(j,0)  =   vt.row(j)  * ((Mt.row(j)).transpose()) ;
            }


            // assign block vector results to final vector (ans) of results
            long  counter = 0;
            for(long j=start_row1; j < start_row1 + num_rows_in_block1; j++){
                 ans(j,0) = ans_tmp(counter,0);
                 var_ans(j,0) = var_ans_tmp(counter,0);
                 counter++;
            }
    
             if (quiet)  Rcpp::Rcout << "block done ... " << std::endl;


      } // end for long



}  //  end if block update


  return Rcpp::List::create(Rcpp::Named("a")=ans,
                            Rcpp::Named("vara") = var_ans);


}






// [[Rcpp::export]]
void createM_ASCII_rcpp(CharacterVector f_name, CharacterVector f_name_ascii, 
                  CharacterVector  type,
                  string AA,
                  string AB, 
                  string BB,
                  double  max_memory_in_Gbytes,  std::vector <long> dims,
                  bool csv, 
                  bool quiet) 
{
  // Rcpp function to create space-removed ASCII file from ASCII and PLINK input files

size_t found;


std::string 
   line; 


ofstream
   fileOUT;

int 
   genoval;

std::string 
     ftype = Rcpp::as<std::string>(type),
     fname = Rcpp::as<std::string>(f_name),
     fnameascii = Rcpp::as<std::string>(f_name_ascii);



//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double 
  memory_needed_in_Gb;
  if (ftype == "PLINK"  ){
  // this is a PLINK ped file. Hence, we need to adjust the dims[1] to get the 
  // size of the genotype file in R land. 
    memory_needed_in_Gb =  (dims[0] *  (dims[1]-6.0)/2.0  *   sizeof(double) )/( (double) 1000000000) ;
  } else {
    // text file
    memory_needed_in_Gb =  (dims[0] *  dims[1] *   sizeof(double) )/( (double) 1000000000) ;
  }
  if ( ftype == "PLINK"  ){
     //------------------------------------
     // convert PLINK ped file into ASCII file with no spaces
     //----------------------------------------------
       CreateASCIInospace_PLINK(fname, fnameascii, dims, quiet);

   }  else {
      //-------------------------------------------
      // convert text file into ASCII file
      //-----------------------------------------
      // Here, we do not need to worry about the amount of memory because 
      // we are processing a line of the file at a time. This is not the case when 
      // creating a ASCII Mt because we have to read in blocks before we can 
      // transpose. 
      if (quiet)
          Rcout << " A text file is being assumed as the input data file type. " << std::endl;
      CreateASCIInospace(fname, fnameascii, dims, AA, AB, BB, csv, quiet);
   }  // end if type == "PLINK" 

//--------------------------------------
// Summary of Genotype File
//--------------------------------------

Rcpp::Rcout <<  "\n\n                    Summary of Marker File  " << std::endl;
Rcpp::Rcout <<  "                   ~~~~~~~~~~~~~~~~~~~~~~~~   " << std::endl;
Rcpp::Rcout <<  " File type:                " << type  << std::endl;
Rcpp::Rcout <<  " File name:                " << fname << std::endl;
Rcpp::Rcout <<  " New ASCII file name:  " << fnameascii  << std::endl;
Rcpp::Rcout <<  " Number of individuals:    "     << dims[0] << std::endl;
if (ftype == "PLINK"  ){
Rcpp::Rcout <<  " Number of loci:           "  << (dims[1] -6)/2.0   << std::endl;
} else {
Rcpp::Rcout <<  " Number of loci:           "  << dims[1] << std::endl;
}
Rcpp::Rcout.precision(2);
Rcpp::Rcout <<  " File size (Gbytes):       "  << memory_needed_in_Gb << std::endl;
Rcpp::Rcout <<  " Available memory (Gbytes):" << max_memory_in_Gbytes  << std::endl;
Rcpp::Rcout << "\n\n" << std::endl;


}











// [[Rcpp::export]]
Eigen::VectorXi  extract_geno_rcpp(CharacterVector f_name_ascii, 
                                   double  max_memory_in_Gbytes, 
                                    long  selected_locus, 
                                    std::vector<long> dims,
                                   Rcpp::NumericVector indxNA)
{
  std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

  long 
     nind;

  nind = dims[0];



  // if (!R_IsNA(indxNA(0)))
  if (indxNA.size() != 0 )
    nind = dims[0] - indxNA.size();


//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
double
  memory_needed_in_Gb =  (dims[0] *  dims[1] *   sizeof(double) )/( (double) 1000000000) ;


Eigen::VectorXi
   column_of_genos(nind);



if(max_memory_in_Gbytes > memory_needed_in_Gb ){
   // reading entire data file into memory
    Eigen::MatrixXd genoMat =  ReadBlock(fnamebin,  0, dims[1], dims[0]);

    // removing rows that correspond to individuals with no 
  // trait data
  if(indxNA.size() != 0){
     for (long ii=0; ii < indxNA.size(); ii++){
           removeRow(genoMat, indxNA(ii) ); } 
  }
   column_of_genos = genoMat.col(selected_locus).cast<int>() ;
   


}  else {
   Rcout << " In else part of if statement " << endl;
    long num_rows_in_block = (max_memory_in_Gbytes  * (double) 1000000000 )/(sizeof(double) * dims[1]);

         long num_blocks = dims[0]/num_rows_in_block;
          if (dims[0] % num_rows_in_block)
                 num_blocks++;


          for(long i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1  ;

              Eigen::MatrixXd    
                genoMat_block1 ( ReadBlock(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;

              // removing rows that correspond to individuals with no 
              // trait data
              long start = start_row1;
              long  finish = start_row1 + num_rows_in_block1;
              // if(!R_IsNA(indxNA(0))){
              if(indxNA.size() != 0 ){
                 for (long ii=0; ii < indxNA.size(); ii++){
                   if(indxNA(ii) >= start & indxNA(ii) <= finish)
                        removeRow(genoMat_block1, indxNA(ii) );
                 }
              }


              // dealing with assigning column_of_genos when some values 
              // may be missing due to having been removed. 
              long colindx = start_row1;
              for(long j=start_row1; j< start_row1+num_rows_in_block1 ; j++){
                 bool found = 0;
                 for(long ii = 0; ii < indxNA.size() ; ii++){
                   if(indxNA[ii] == j){
                          found=1;
                   }
                 }
                 if (!found){ 
                    // j not in indxNA
                   column_of_genos(colindx) = genoMat_block1.col(selected_locus)(j-start_row1);
                   colindx++;
                }

              } // end for j

          } // end for  i


} // end if max_memory

return(column_of_genos);

}






// [[Rcpp::export]]
Eigen::MatrixXd  calculateMMt_rcpp(CharacterVector f_name_ascii, 
                                   double  max_memory_in_Gbytes, int num_cores,
                                   Rcpp::NumericVector  selected_loci , std::vector<long> dims, 
                                   bool quiet) 
{
// set multiple cores
Eigen::initParallel();
omp_set_num_threads(num_cores);
Eigen::setNbThreads(num_cores);
Rcpp::Rcout << " Number of cores being used for calculation is .. "  << Eigen::nbThreads() << endl;


std::string 
   line; 


ofstream
   fileOUT;

int 
   genoval;

std::string 
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

// gpu will only work with double precision matrices in Eigen. 
// Had to change code to be double precision. 
MatrixXd
    MMt(MatrixXd(dims[0], dims[0]).setZero());

//MatrixXi 
//    genoMat,
//    MMt(MatrixXi(dims[0], dims[0]).setZero());





//-----------------------------------
// Calculate amount of memory needed
//-----------------------------------
// Memory required for 
// MMt   dims[0] * dims[0] *  sizeof(double) 
// genoMat  dims[0] * dims[1] * sizeof(double) 
// genoMat transpose dims[0] * dims[1] * sizeof(double) 
// 
// Block update
//
// MMt   dims[0] * dims[0] *  sizeof(double) 
// genoMat  num_rows_in_block * dims[1] * sizeof(double) 
// genoMat transpose num_rows_in_block * dims[1] * sizeof(double) 

double 
  memory_needed_in_Gb =  (dims[0]*dims[0]* sizeof(double)  + 2*(dims[0] *  dims[1] *   sizeof(double) ))/( (double) 1000000000) ;




//-------------------------
// Perform MMt calculation
//-------------------------
if(max_memory_in_Gbytes > memory_needed_in_Gb ){
   // reading entire data file into memory
    Eigen::MatrixXd genoMat = ReadBlock(fnamebin,  0, dims[1], dims[0] );
   if(!R_IsNA(selected_loci(0))){
     // setting columns to 0
     for(long ii=0; ii < selected_loci.size() ; ii++) 
       genoMat.col(selected_loci(ii)).setZero(); 
   }


   MMt.noalias() = genoMat * genoMat.transpose(); 



} else {
    // based on user defined memory. Doing MMt via blockwise multiplication
    // long num_rows_in_block = (max_memory_in_Gbytes  * (double) 1000000000 )/(sizeof(double)  * dims[1]);
 //   long num_rows_in_block = (max_memory_in_Gbytes  * 1000000000 - dims[0] * dims[0] * sizeof(double) )/( 2* sizeof(double)  * dims[1]);
   // its 2.2 instead of 2 and 4.84 instead of 4 to give us a 10% memory buffuer
    double part1 = -2.2 *  dims[1];
    double part2 = 4.84*dims[1] * dims[1] +  4 * max_memory_in_Gbytes  * 1000000000.0/sizeof(double);
    part2 = sqrt(part2);
    long num_rows_in_block = (part1 + part2)/2.0;
    Rcout << "number of rows in block is " << num_rows_in_block << endl;  

           // blockwise multiplication

          // find out the number of blocks needed
          long num_blocks = dims[0]/num_rows_in_block;



          if (dims[0] % num_rows_in_block)
                 num_blocks++;


          for(long i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1;
            //  Rcpp::Rcout << num_rows_in_block1 << " num rows in block 1 " << std::endl;


               Rcout << " Reading Block " << i << "  Data .... ------------- " << endl;
               Eigen::MatrixXd    
                    genoMat_block1 ( ReadBlock(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;
               Rcout << " Finished reading Block Data .... --------------- " << endl;

              Rcout << " taking MMtsub(MatrixXd(num_rows_in_block1, num_rows_in_block1).setZero()) " << endl;
              Eigen::MatrixXd    
                   MMtsub(MatrixXd(num_rows_in_block1, num_rows_in_block1).setZero());

             if(!R_IsNA(selected_loci(0) )){
             // setting columns to 0
             for(long ii=0; ii < selected_loci.size() ; ii++)
                genoMat_block1.col(selected_loci(ii)).setZero();
             }
             // Rcpp::Rcout << "  Block 1  "  << std::endl;
             // Rcpp::Rcout << genoMat_block1.rows() << std::endl;
             // Rcpp::Rcout << genoMat_block1.cols() << std::endl;

              Rcout << "---------------------GPU--C++ -------------------------   " << endl;
              Rcout << " MMtsub = genoMat_block1 * genoMat_block1.transpose(); " << endl;
              Rcout << "------------------------------------------------------ " << endl;
              MMtsub.noalias() = genoMat_block1 * genoMat_block1.transpose(); 

              //          i            j            num rows               num   cols
              MMt.block(start_row1, start_row1, num_rows_in_block1, num_rows_in_block1) = MMtsub;


              for(long j=i+1;j<num_blocks; j++){
                   long start_row2 = j * num_rows_in_block;
                   long num_rows_in_block2 = num_rows_in_block;
                   if ((start_row2 + num_rows_in_block2) > dims[0])
                          num_rows_in_block2 = dims[0] - start_row2;
                    Rcout << " Reading genoMat_block " << j << "  data ... "   << endl;
                    Eigen::MatrixXd    
                       genoMat_block2 ( ReadBlock(fnamebin,  start_row2, dims[1], num_rows_in_block2)) ;




                   Eigen::MatrixXd    MMtsub(MatrixXd(num_rows_in_block1, num_rows_in_block2).setZero());

                  if(!R_IsNA(selected_loci(0) )){
                   // setting columns to 0
                   for(long jj=0; jj < selected_loci.size() ; jj++)
                      genoMat_block2.col(selected_loci(jj)).setZero();
                   }
            //  Rcpp::Rcout << " Block 2 " << std::endl;
            //  Rcpp::Rcout << genoMat_block2.rows() << std::endl;
            //  Rcpp::Rcout << genoMat_block2.cols() << std::endl;
              Rcout << "---------------------GPU--C++ -------------------------   " << endl;
              Rcout << " MMtsub = genoMat_block1 * genoMat_block2.transpose(); " << endl;
              Rcout << "------------------------------------------------------ " << endl;
                   MMtsub.noalias() = genoMat_block1 * genoMat_block2.transpose(); 
                   //          i,        j,     num rows,              num cols
                   MMt.block(start_row1, start_row2, num_rows_in_block1, num_rows_in_block2) = MMtsub;
                   // and its symmetric block
                   MMt.block(start_row2, start_row1, num_rows_in_block2, num_rows_in_block1) = MMtsub.transpose();


            }  // end for int j





          } // end for int


 // }  // end inner if else

}  // end outer if else



  return MMt;

}




