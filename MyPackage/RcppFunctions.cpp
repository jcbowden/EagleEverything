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
//          5. read blocks of binary to form sub matrices of M.
//          6. Perform M %*% t(M) as a block multiplication. 
//

// This was causing issues when building on clean Linux system
// creates reliance on mkl.h 
#define EIGEN_USE_BLAS

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <R.h>

// added by Ryan
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <time.h>
// end added by Ryan

#include <omp.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <bitset>
#include <string>
#include <fcntl.h>
#include <stdlib.h>
// #include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>

//    #include <magma.h>



using namespace std;
using namespace Rcpp;
using namespace RcppEigen;

using Eigen::MatrixXi;
using Eigen::MatrixXd;  
using Eigen::Lower;
using Eigen::Map;   // maps rather than copies

#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif


// const size_t bits_in_double = std::numeric_limits<long double>::digits;
const size_t bits_in_ulong = std::numeric_limits<unsigned long int>::digits;
const size_t bits_in_int = std::numeric_limits<int>::digits;



// [[Rcpp::export]]
std::vector <long>    ReshapeM_rcpp( CharacterVector  fnameM, 
               CharacterVector  fnameMt, 
               std::vector <long> indxNA, 
               std::vector <long> dims){


   // note indxNA starts from 0


  std::vector <long> 
        newdims(2,0);

 std::ostringstream
      os;



   std::string
       line, 
       FnameM = Rcpp::as<std::string>(fnameM),
       FnameMt = Rcpp::as<std::string>(fnameMt); 



   //-------------------------------------------
   // converting M.ascii to reshaped M.asciitmp
   //-------------------------------------------
   // open file and check for its existence. 
   std::ifstream fileIN(FnameM.c_str());
   if(!fileIN.good()) {
        os << "\n\nERROR: Could not open  " << FnameM << "\n\n" << std::endl;
        Rcpp::stop(os.str() );
   }


   // change name for new no-space ASCII file with rows matching indxNA removed
   FnameM.append("tmp");


   // open ascii file that is to hold no-spaces genotype data
   std::ofstream fileOUT(FnameM.c_str(), std::ios::out );

 long rownum=0;
 bool writeline;
 while(fileIN.good()){
      while(getline(fileIN, line)){
          writeline = true;
          for(unsigned long ii=0; ii<indxNA.size(); ii++){
            if(indxNA[ii] == rownum)
                writeline = false;
          }
          if(writeline){
              fileOUT << line << endl;
              newdims[0]++;
          }
          rownum++;

      newdims[1] = line.length();

      }  // end inner while

 }  // end outer while(fileIN

 fileIN.close();
 fileOUT.close();

  //-------------------------------------------
   // converting Mt.ascii to reshaped Mt.asciitmp
   //-------------------------------------------

   // open file and check for its existence. 
   std::ifstream fileINt(FnameMt.c_str());
   if(!fileINt.good()) {
        os << "\n\nERROR: Could not open  " << FnameMt << "\n\n" << std::endl;
        Rcpp::stop(os.str() );
   }


   // change name for new no-space ASCII file with rows matching indxNA removed
   FnameMt.append("tmp");

 
   // open ascii file that is to hold no-spaces genotype data
   std::ofstream fileOUTt(FnameMt.c_str(), std::ios::out );

  while(fileINt.good()){
      while(getline(fileINt, line)){
          // removing columns
          // this is okay since indxNA is in decreasing size
          for(unsigned long ii=0; ii<indxNA.size(); ii++){
            line.erase ( indxNA[ii], 1 );
          }  // end for a
          fileOUTt << line << endl;
      }  // end inner while
 }  // end outer while(fileIN

 fileINt.close();
 fileOUTt.close();
  return newdims;

}  //end ReshapeM





// ---- code developed by Ryan
//
//char* mapFileFromDiscBlocked(const char * file_name, unsigned long long &sizeUsed, unsigned long long &sizeActual, unsigned long long blockSize, unsigned long long offset) {
//	// If set to true debugging information is displayed
//	// Otherwise these messages are suppressed
//	bool debugMsgs = true;
//
//	// Read file size and system page file size
//	int pagesize = getpagesize();
//
//	// Open file descriptor for given file name
//	// Read-only permission
//	int fd;
//
//	// Open file and get file descriptor
//	try {
//		fd = open (file_name, O_RDONLY);
//		// File descriptor returned as -1 if file cannot be opened
//		if ( -1 == fd) throw "File could not be opened for reading";
//	}
//	catch (const char* msg) {
//                Rcpp::stop(msg);
//                return "";
//
//	}
//
//	// Round up file size to next multiple
//	// of system page size;
//	sizeActual = blockSize;
//	sizeUsed = sizeActual + (pagesize - (sizeActual % pagesize));
//	// if (debugMsgs) printf("Memory used to map file: %2.5f Mbytes\n", sizeUsed/1024/1024);
//
//	// Memory mapped data file
//	char *fileMemMap;
//
//	// Map file to memory using mmap() system function
//	// File may now be treated as a character array in memory
//	// Syntax:
//	// 		void *mmap(void *addr, size_t length,
//	//			  int prot, int flags, int fd, off_t offset);
//	// Permissions:
//	//		PROT_READ: Pages in memory only allow read-only operations
//	//		MAP_PRIVATE: Changes not visible to other processes
//	//					 Underlying file is no altered
//
//	// cout << "Page Size " << pagesize << endl;
//	// cout << "Size Used " << sizeUsed << endl;
//	// cout << "Size Actual " << sizeActual << endl;
//
//	fileMemMap = (char *) mmap (0, sizeUsed, PROT_READ, MAP_PRIVATE, fd, offset);
//
//	//if (madvise(fileMemMap, sizeUsed-1024, MADV_WILLNEED | MADV_SEQUENTIAL) == -1) {
//	//	cout << "madvise error" << endl;
//	//	return NULL;
//	//}
//
//	// Close the file descriptor
//	// Frees resources associated with the file descriptor
//	close(fd);
//
//	return fileMemMap;
//}
//
//char* mapFileFromDisc(const char * file_name, unsigned long &sizeUsed, unsigned long &sizeActual, Rcpp::Function message) {
//	// If set to true debugging information is displayed
//	// Otherwise these messages are suppressed
//	bool debugMsgs = true;
//
//	// Read file size and system page file size
//	int pagesize = getpagesize();
//
//	// Stores information about file
//	struct stat s;
//
//	// Open file descriptor for given file name
//	// Read-only permission
//	int fd;
//
//	// Open file and get file descriptor
//	try {
//		fd = open (file_name, O_RDONLY);
//		// File descriptor returned as -1 if file cannot be opened
//		if ( -1 == fd) throw "File could not be opened for reading";
//	}
//	catch (const char* msg) {
//		message( "ERROR occurred in mapFileFromDisc " , msg );
//		return NULL;
//	}
//
//	// Get the file size on disc
//	int status = fstat (fd, &s);
//
//	// Round up file size to next multiple
//	// of system page size;
//	sizeActual = s.st_size;
//	sizeUsed = sizeActual + (pagesize - (sizeActual % pagesize));
//	if (debugMsgs) message("Memory used to map file: %d Mbytes\n", sizeUsed/1024/1024);
//
//	// Memory mapped data file
//	char *fileMemMap;
//
//	// Map file to memory using mmap() system function
//	// File may now be treated as a character array in memory
//	// Syntax:
//	// 		void *mmap(void *addr, size_t length,
//	//			  int prot, int flags, int fd, off_t offset);
//	// Permissions:
//	//		PROT_READ: Pages in memory only allow read-only operations
//	//		MAP_PRIVATE: Changes not visible to other processes
//	//					 Underlying file is no altered
//
//	fileMemMap = (char *) mmap (0, sizeUsed, PROT_READ, MAP_PRIVATE, fd, 0);
//
//	if (madvise(fileMemMap, sizeUsed, MADV_WILLNEED | MADV_SEQUENTIAL) == -1) {
//		message( "madvise error" );
//		return NULL;
//	}
//
//	// Close the file descriptor
//	// Frees resources associated with the file descriptor
//	close(fd);
//
//	return fileMemMap;
//}
//
//
//
//
//
//bool  CreateASCIInospaceFast(std::string fname, std::string asciifname, std::vector<long> dims,
//		std::string  AA,
//		std::string AB,
//		std::string BB,
//		bool quiet, 
//                Rcpp::Function message, 
//                std::string missing)
//{
//
//	// Used to store size of memory used for mapping file
//	// Size is rounded up to memory page size, hence two variables
//	// SizeUsed is used to ummap the mapped memory
//	// SizeActual is used to determine the number of characters in the file
//	unsigned long sizeUsed = 0;
//	unsigned long sizeActual = 0;
//	// Map file from hard-disk to memory
//	char* dataFile = mapFileFromDisc(fname.c_str(), sizeUsed, sizeActual, message);
//
//	// Check that an error did not occur in the file mapping process
//	if ( dataFile == NULL )
//	{
//	        message("Error mapping file.");
//		return false;
//	}
//
//	// Array used for buffering file output
//	const long long bufferSize = 8*1024*1024;
//
//
//
//	// char outputBuffer[bufferSize];
////        char* outputBuffer = NULL;
////        outputBuffer = new char[bufferSize];
//
//std::vector<char> outputBuffer (bufferSize);  // allocated on stack, with a data buffer that is probably on the heap
// 
//
//	int inc = 0;
//
//	// Output file
//	FILE *outputFile;
//	outputFile = fopen (asciifname.c_str(), "wb");
//
//	// Assuming BA is the reverse of AB
//	std::string BA = std::string ( AB.rbegin(), AB.rend() );
//	// Find a more elegant way to do this
//	const short AALen = AA.length();
//	const short ABLen = AB.length();
//	const short BBLen = BB.length();
//	int maxLenGen = 0;
//	if ( AALen > ABLen )
//	{
//		maxLenGen = AALen;
//	}else
//	{
//		if ( BBLen > ABLen )
//		{
//			maxLenGen = BBLen;
//		} else {
//			maxLenGen = ABLen;
//		}
//	}
//
//	int slidingInc = 0;
//	 char windowBuffer[maxLenGen];
//
//
//        if (!quiet){
//              message("");
//              message(" Reading text File  ");
//              message("");
//              message(" Loading file ");
//       }
//
//
//
//	int latch = 0;
//	// Loop through file in memory
//	for (unsigned long long i = 0; i <= sizeActual; i++ ) {
//
//		if ( dataFile[i] != ' ' && dataFile[i] != '\n' && i != sizeActual)
//		{
//			latch = 0;
//			windowBuffer[slidingInc] = dataFile[i];
//			slidingInc++;
//		}else {
//			if ( latch == 0 )
//			{
//				// Output magic
//				windowBuffer[slidingInc] = '\0';
//				// cout << "(" << windowBuffer << ")" << endl;
//
//				// Comparison magic
//				if ( AA == windowBuffer ) {
//					// cout << "AA Found" << endl;
//					outputBuffer[inc] = '0';
//					inc++;
//				} else if ( AB == windowBuffer ) {
//					// cout << "AB Found" << endl;
//					outputBuffer[inc] = '1';
//					inc++;
//				} else if ( BA == windowBuffer ) {
//					// cout << "BA Found" << endl;
//					outputBuffer[inc] = '1';
//					inc++;
//				} else if ( BB == windowBuffer ) {
//					// cout << "BB Found" << endl;
//					outputBuffer[inc] = '2';
//					inc++;
//                               } else if ( missing == windowBuffer ) {
//                                   // setting any missing values to het code
//                                        outputBuffer[inc] = '1';
//                                        inc++;
//				} else {
//                                    std::string str(windowBuffer);
//                                    if (AB=="NA"){
//                                       message( "Marker file contains marker genotypes that are different to AA=" , AA , " BB=" , BB);
//                                       message(" For example  " ,  str);
//                                       message(" ReadMarker has terminated with errors");
//                                       return false;
//                                   } else {
//                                       message( "Marker file contains marker genotypes that are different to AA=" , AA , " AB=" , AB , " BB=" , BB);
//                                       message( "For example , " , str );
//                                       message( "ReadMarker has terminated with errors");
//                                       return false;
//                                  }
//
//				}
//
//				// Reset sliding index
//				slidingInc = 0;
//				latch = 1;
//			}
//			if ( dataFile[i] == '\n' )
//			{
//				outputBuffer[inc] = '\n';
//				inc++;
//			}
//		}
//
//		// If output buffer is full
//		if ( inc >= bufferSize)
//		{
//			// Write buffer to file on disk
//                        for(long ii=0; ii < outputBuffer.size(); ii++)
//	                   fwrite (&outputBuffer[ii] , sizeof(char), 1 , outputFile);
//			//cout << inc << endl;
//			// Reset buffer index
//			inc = 0;
//		}
//	}
//
//	// Writing remaining data in output buffer to output file
//
//        for(long ii=0; ii < inc; ii++)
//          fwrite (&outputBuffer[ii] , sizeof(char), 1 , outputFile);
//
//	// No longer need to use memory mapped file, release it
//	munmap(dataFile, sizeUsed);
//
//
//
//  // write out a few lines of the file if quiet
////  if(quiet > 0){
//     // open ascii  file
//     std::string line, tmp;
//     std::ifstream fileIN(fname.c_str());
//     long counter = 0;
//     int nrowsp =  5;
//     int ncolsp = 12;
//     if(dims[0] < 5)
//         nrowsp = dims[0];
//     if(dims[1] < 12)
//          ncolsp = dims[1];
//
//     message(" First ", nrowsp, " lines and ", ncolsp, " columns of the marker file. ");
//     std::string rowline;
//     while(getline(fileIN, line ) && counter < nrowsp)
//     {
//       std::ostringstream oss;
//       std::istringstream streamA(line);
//       for(int i=0; i < ncolsp ; i++){
//           streamA >> tmp;
//           oss << tmp << " " ;
//        }
//        std::string rowline = oss.str();
//        message(rowline);
//        counter++;
//      }  // end  while(getline(fileIN, line ))
////  } // end if(quiet)
//
//
//
//
//
//
//	// Close output file
//	fclose (outputFile);
//	return true;
//}
//
//
//
//
//
//
//// ----- end of code by Ryan
//
//
//
//
//
//
//
//std::vector<std::string> split(const char *str, char c = ' ')
//{
//    std::vector<std::string> result;
//
//    do
//    {
//        const char *begin = str;
//
//        while(*str != c && *str)
//            str++;
//
//        result.push_back(std::string(begin, str));
//    } while (0 != *str++);
//
//    return result;
//}
//
//






//get number of rows and columns in marker file
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// [[Rcpp::export]]
std::vector<long>   getRowColumn(std::string fname) 
{
  // Purpose:  to open the marker file where the marker data are kept.
  //           An error will be produced if the file cannot be found.
  //           I am assuming no row or column names
// int 
//   genoval;

 std::string
   line;

 std::ostringstream 
      os;


 std::vector<long> dimen(2,0)  ;  // dim[0] row number
                               // dim[1] col number 

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
 fileIN.seekg(0, std::ios::beg);

 getline(fileIN, line );
 std::istringstream streamA(line);

  std::string 
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
                         bool quiet, Rcpp::Function message)
{
//long 
//   colindx = 0; 

int
  n_of_cols_in_geno = (dims[1] -6)/2.0;






// char alleles [ 2 ][ n_of_cols_in_geno ];  // holds alleles  
//char* alleles0 = NULL;
//      alleles0 = new char[ n_of_cols_in_geno ];
//char* alleles1 = NULL;
//      alleles1 = new char[ n_of_cols_in_geno ];
std::vector<char> alleles0( n_of_cols_in_geno );
std::vector<char> alleles1( n_of_cols_in_geno );


std::vector<char>
     rowvec( dims[1] - 6 );  // holds allelic information from PLINK file


std::string
   tmp,
   token,
   line;

//char 
//   sep = ' ';


 std::ostringstream 
      os;




//------------------------------

// open PLINK ped  file
std::ifstream fileIN(fname.c_str());
if(!fileIN.good()) {
  message("\nERROR: PLINK ped file could not be opened with filename  ",   fname );
  os << "\n\nERROR: ReadMarkerData has terminated with errors.  " << fname << "\n\n" << std::endl;
  Rcpp::stop(os.str() );
}

// open ascii file that is to hold no-spaces genotype data
std::ofstream fileOUT(asciifname.c_str(), std::ios::out );
long  counter = 0;

// initializing input line 
std::string rowinfile(n_of_cols_in_geno, '0'); // s == "000000"

int printOnlyOnce = 0;  // flag for printing warning message about missing data


while(getline(fileIN, line ))
{


  std::istringstream streamLine(line);

 // check number of columns for each line
 long number_of_columns = 0;
 std::string rowinfile(n_of_cols_in_geno, '0'); // s == "000000"
 
 if (!quiet ){
    while(streamLine >> tmp)
        number_of_columns ++;

     if (number_of_columns != dims[1] ){
         message("\n");
         message( "Error:  PLINK file contains an unequal number of columns per row.  " );
         message( "        The error has occurred at row " ,  counter+1 , " which contains " ,  number_of_columns ,  " but " );
         message( "        it should contain " , dims[1] , " columns of data. " );
         message("\n");
          os << " ReadMarkerData has terminated with errors\n" << std::endl;
         Rcpp::stop(os.str() );
       }  // end  if (number_of_columns != dims[1] )
   } // end  if (quiet )

   std::istringstream streamA(line);
   // tokenized row and placed it in std::vector rowvec
   for(int i=0; i <= 5; i++)
     streamA >> tmp;
   for(long i=6; i < dims[1] ; i++){
            streamA >> rowvec[i-6];
   }  // end  for(long i=0; i < dims[1] ; i++)


   // initialize alleles structure to first row of PLINK info
   if (counter == 0) {
         for(long i=0; i < n_of_cols_in_geno ; i++){
            if( rowvec[ (2*i ) ] == '0' ||  rowvec[ (2*i + 1) ] == '0' || rowvec[ (2*i ) ] == '-' ||  rowvec[ (2*i + 1) ] == '-'){
               // missing allele
               alleles0[ i ] = 'I';
               alleles1[ i ] = 'I';
            } else {
               alleles0[ i ] =  rowvec[ (2*i ) ];
               alleles1[ i ] =  rowvec[ (2*i + 1) ];
            } //end  if (rowvec 
         }
   }

   // turn allelic info from PLINK into genotype 0,1,2 data
     // also do some checks for more than 2 alleles, and 0 and - for missing data
     for(long i=0; i < n_of_cols_in_geno; i++){
        // Checking for missing allelic information in PLINK file
        if( rowvec[ (2*i ) ] == '0' ||  rowvec[ (2*i + 1) ] == '0' || rowvec[ (2*i ) ] == '-' ||  rowvec[ (2*i + 1) ] == '-'){
           if (printOnlyOnce == 0){
                 message("\n");
                 message(" Warning:  PLINK file contains missing alleles (i.e. 0 or - ) " );
                 message("           These missing genotypes should be imputed before running Eagle." );
                 message("           As an approximation, AMpus has set these missing genotypes to heterozygotes. " );
                 message("           Since Eagle assumes an additive model, heterozygote genotypes do not contribute to the estimation of " );
                 message("           the additive effects.  " );
                 message("\n");
                 printOnlyOnce = 1;
            } // if printOnlyOnce`
            rowvec[ (2*i) ] = 'I';      // impute
            rowvec[ (2*i + 1) ] = 'I';  // impute
        }


        // Check if allele has been seen before in allele file. 
        // If so, make sure alleles doesn't already  contain two alleles - otherwise generate error message
        for(int j = 1; j >= 0; --j){ // looping over the two alleles with indexes 0 and 1
           if (rowvec[ (2*i + j) ] != alleles0[ i ] && rowvec[ (2*i + j) ] != alleles1[ i ]){
              // situation 1: rowvec contains missing values ie 'I' then do nothing
              if (rowvec[ (2*i + j) ] == 'I' ){
                 // do nothing here

              } else {
                // situation 2: alleles contain missing values I
                if (alleles0[i] == 'I'){
                     alleles0[i] = rowvec[ (2*i + j) ];
                } else {
                      if (alleles1[i] == 'I'){
                           alleles1[i] = rowvec[ (2*i + j) ];
                       } else {
                         if (alleles0[ i ] == alleles1[ i ] ){
                             // this is okay. alleles only contains a single allele at the moment. Re-initialise alleles
                             alleles1[ i ] = rowvec[ (2*i + j) ];
                           } else {
                              // Error - we have more than two alleles segregating at a locus
                            message("\n");
                            message("Error:  PLINK file cannot contain more than two alleles at a locus.");
                            message("        The error has occurred at snp locus " , i  + 1 , " for individual " , counter+1 );
                            message("\n");
                            os << " ReadMarkerData has terminated with errors\n" << std::endl;
                             Rcpp::stop(os.str() );
                          } // end inner if else

                       } // end if (alleles1[i] == 'I')
                } // end  if (alleles0[i] == 'I')

              } // end if (rowvec[ (2*i + j) ] == 'I' )


           }  // end if (rowvec[ (2*i + j) ] != alleles0[ i ] && rowvec[ (2*i + j) ] != alleles1[ i ])




    // set rowinfile
    if (rowvec[ (2*i) ] == 'I' || rowvec[ (2*i + 1) ] == 'I' ){
        rowinfile[i] = '1' ; // AB geno no additive effect
    } else {
        if (rowvec[ (2*i + 1) ] !=   rowvec[ (2*i) ] ){
          rowinfile[i] = '1' ;  // AB
        } else {
          if (rowvec[ (2*i ) ] == alleles0[ i ] ){  // matches first allele
               rowinfile[i] = '0';  // AA
          }  else {
               rowinfile[i] = '2';  // BB
          }
        }  // end outer if else rowvec
    } // end  if (rowvec[ (2*i) ] == 'I' || rowvec[ (2*i + 1) ] == 'I' )
 } // end  for(long i=0; i < n_of_cols_in_geno; i++)




  }  // end for(long i=0; i< n_of_cols_in_geno ; i++)

  fileOUT << rowinfile;
  fileOUT << "\n";
  counter++;


  }  // end while(getline(fileIN, line ))



// write out a few lines of the file if quiet
// open PLINK ped  file
std::ifstream fileIN_backtobeginning(fname.c_str());
counter = 0;


int nrowsp =  5;
int ncolsp = 24;
if(dims[0] < 5)
nrowsp = dims[0];
if(dims[1] < 25)
ncolsp = dims[1];

message(" First ", nrowsp, " lines and ", ncolsp, " columns of the PLINK ped file. ");
 

std::string rowline;

while(getline(fileIN_backtobeginning, line ) && counter < nrowsp)
{
       std::ostringstream oss;
       std::istringstream streamB(line);
       for(int i=0; i < ncolsp ; i++){
           streamB >> tmp;
           oss << tmp << " " ;
        }
        std::string rowline = oss.str();
        message(rowline);
        counter++;
}  // end  while(getline(fileIN, line ))




// close files
fileIN.close();
fileOUT.close();


}











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

//long 
//   colindx = 0;




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
// std::string rowinfile(dims[1], '0'); // s == "000000"

 // std::string* rowinfile = NULL;
 // rowinfile = new std::string  [ dims[1] ];

std::vector<char> rowinfile( dims[1] );


 for(long i=0; i < dims[1]; i++)
    rowinfile[i] = '0';



while(getline(fileIN, line ))
{


   Rcout << "\r" << 100.0*counter/dims[0] << "% read of text file.       " << flush;
 // message(100.0*counter/dims[0] ,  "% read of text file.       ");

 // Here, BB is coded into 2 
 //       AB is coded into 1, 
 //       AA is coded into 0. 
  std::istringstream streamA(line);
  long i=0;
  number_of_columns = 0;
  while(streamA >> token)
  {
 //    if(quiet )
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
              message( "Marker file contains marker genotypes that are different to AA=" , AA , " BB=" , BB);
              message(" For example , " , token );
              message(" ReadMarker has terminated with errors");
              return false;
          } else {
              message( "Marker file contains marker genotypes that are different to AA=" , AA , " AB=" , AB , " BB=" , BB);
              message( "For example , " , token );
              message( "ReadMarker has terminated with errors");
              return false;
         }
       }  //end if else 
       i++;
  } // end whle streamA

  if (!quiet){
        message(" Number of columns in line " , counter+1 , " is " , number_of_columns);

        if (number_of_columns != dims[1] ){
             message("\n");
             message("Error:  Marker text file contains an unequal number of columns per row.  ");
             message("        The error has occurred at row " , counter+1 , " which contains " , number_of_columns , " but ");
             message("        it should contain " , dims[1] , " columns of data. ");
             message("\n");
             message(" ReadMarkerData has terminated with errors");
             return false;
       }  // end if number_of_columns
  } // end if quiet
  for(long ii=0; ii< number_of_columns; ii++){
     fileOUT << rowinfile[ii];
  }
  fileOUT << "\n";
  counter++;

 }  // end while getline









  // write out a few lines of the file if quiet
//  if(quiet ){
     // open ascii  file
//     line;
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
//  } // end if(quiet)

// close files
fileIN.close();
 fileOUT.close();
// fclose(pfileOUT);

  return true;
}




//
//// Ryan's ReadBlock code which uses mmap() system call
//// Not fully tested, use with caution
//Eigen::MatrixXd  ReadBlockFast(std::string asciifname,
//		long start_row,
//		long numcols,
//		long numrows_in_block) {
//
//	// Start settable parameters
//	// const unsigned long long maxMemory = 16ull*1024*1024*1024; // Mbytes
//	// End settable parameters
//
//	long long pagesize = getpagesize();
//	unsigned long long sizeUsed = 0;
//	unsigned long long sizeActual = 0;
//
//	// Used for offseting to a newline from a loaded block
//	unsigned long long realPos = (start_row)*(numcols+1);
//	unsigned long long allignedPos = realPos - (realPos % pagesize);
//	unsigned long long offsetCol = (realPos - allignedPos) % numcols;
//	unsigned long long offsetRow = floor((realPos - allignedPos) / numcols);
//
//	// Debugging information
//	/*cout << "NUM COLS: " << numcols << endl;
//	/cout << "NUM ROWS: " << numrows_in_block << endl;
//	cout << "REAL: " << realPos << endl;
//	cout << "ALLIGNED: " << allignedPos << endl;
//	cout << "Offset Col: " << offsetCol << endl;
//	cout << "Offset Row: " << offsetRow << endl;
//	cout << "Eigen Dimensions: (" << numrows_in_block << ", " << numcols << ")" << endl;
//	cout << "Starting Point: " << offsetRow*(numcols) + offsetCol << endl;
//	cout << "To Read: " << numrows_in_block * (numcols+1) << endl;
//	*/
//	// Eigen matrix to store block of data in
//	Eigen::MatrixXd M(numrows_in_block, numcols) ;
//
//	// Offset point to start reading mapped file from
//	unsigned long startingPoint = offsetRow*(numcols) + offsetCol;
//
//	// Memory required to read block of data
//	unsigned long long requiredMemory = startingPoint + ( numrows_in_block*(numcols+1)) * 8; // bytes
//
//
//	// Memory usage check
//// 	if (requiredMemory <= maxMemory){ 
//
//		// Memory map file
//		char* dataFile = mapFileFromDiscBlocked(asciifname.c_str(), sizeUsed, sizeActual, requiredMemory, allignedPos);
//
//		// Used to load through the mapped file
//		unsigned long long rowInc = 0;
//		unsigned long long colInc = 0;
//
//		// Loop through file data, reads data in as column major (Default for Eigen)
//		for (unsigned long long int i = 0; i < numrows_in_block * (numcols+1); i++) {
//			// If newline character is encountered
//			if ( dataFile[startingPoint+i] == '\n') {
//				colInc++;
//				rowInc = 0;
//			} else {
//				// Read value from data file and convert to number format and subtract one
//				// 0 -> -1, 1 -> 0, 2 -> 1 as per original ReadBlock code
//				signed int value = (dataFile[startingPoint+i]  - '0') - 1;
//
//				// Store value in Eigen matrix
//                                //  Rcout << "colInc= " << colInc << " rowInc " << rowInc << " value = " << value << endl;
//				 M(colInc,rowInc) = value;
//
//				// Used for debugging purposes
//				// cout << "Row: " << rowInc << " :: " << "Col: " << colInc << " :: " << value << endl;
//
//				rowInc++;
//			}
//               } 
//		// Done using memory mapped file, release it
//		munmap(dataFile, sizeUsed);
//
////	} else {
////   cout << "ERROR: Insufficent memory allocated to reading block of data." << endl;
////		cout << "Required: " << requiredMemory/1024/1024/1024 << endl;
////		exit(1);
//// 	}
//
//	return M;
//}
//
//
//
//

Eigen::MatrixXd  ReadBlock(std::string asciifname, 
                           long start_row,
                           long numcols,
                           long numrows_in_block)

{
 // reads in data from ASCII file 
 // to form M Eign double matrix 
std::ostringstream
      os;
std::string
   line;


long 
  // coli = 0, 
  rowi = 0;

//long 
//    igeno;

//double 
//     geno;

Eigen::MatrixXd
      M(numrows_in_block, numcols) ;


// Open no-space ASCII file
   std::ifstream fileIN(asciifname.c_str(), std::ios::in );

    if(!fileIN.good()) {
      os << "ERROR: Could not open  " << asciifname << std::endl;
      Rcpp::stop(os.str() );
     }

   for(long rr=0; rr < (start_row + numrows_in_block) ; rr++){
      // read a line of data from ASCII file
      getline(fileIN, line);
      if(rr >= start_row){
          std::istringstream streamA(line);
          for(long ii=0; ii < numcols  ; ii++){
            int tmp  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
            M(rowi, ii) = (double) tmp - 1;   // converting data to -1, 0, 1 
          }
          rowi++;
      } // end if rr
   } // end for(rr



// Close the ascii file
   fileIN.close();


 return M;

}








// [[Rcpp::export]]
void  createMt_ASCII_rcpp(CharacterVector f_name, CharacterVector f_name_ascii, 
                              double  max_memory_in_Gbytes,  std::vector <long> dims,
                              bool  quiet, Rcpp::Function message )
{

// read data from M.ascii that has already been created and transpose this file

std::string
   token, 
   line;

std::ostringstream
      os;

long
 n_of_cols_to_be_read;

//int
//   genoval;

std::string
     fname = Rcpp::as<std::string>(f_name),
     fnameascii = Rcpp::as<std::string>(f_name_ascii);


//char 
//   sep = ' ';


double 
   max_mem_in_bytes  =  max_memory_in_Gbytes * 1000000000;


// Calculate number of columns that can be read in as a block with XGb of
// memory. 
   // how much memory will be needed to  store M, takes it transpose M.transpose, and 
   // store its answer in Mt (+ .5 for a buffer)
double mem_bytes = 3.5 * dims[0] * dims[1] * (bits_in_int/8);  // assumes a 64 bit system


 // open  files
std::ofstream fileOUT(fnameascii.c_str(), std::ios::out);
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
    fileIN.seekg(0, std::ios::beg);

   // read values into matrix
   long rowi=0;
   while(getline(fileIN, line ))
   {
       // read a line of data from ASCII file
       for(long coli=0; coli < dims[1]; coli++){
           M(rowi, coli)  = line[coli] - '0'; // trick to removes ASCII character offset for numbers
       }
       rowi++;
   }  // end while getline

  // take transpose of matrix M
  MatrixXi Mt = M.transpose();
  
  // write out contents fo Mt to file (no spaces)a
 // std::string rowinfile(Mt.cols(), '0');  // initialising string with 0's
 // std::string* rowinfile = NULL;
 // rowinfile = new string  [ Mt.cols()  ];

 std::vector<char> rowinfile( Mt.cols() );
 for(long i=0; i < Mt.cols() ; i++)
    rowinfile[i] = '0';



  for(long rowi=0; rowi<Mt.rows(); rowi++){
     for(long coli=0; coli<Mt.cols(); coli++){
         // fileOUT << Mt(rowi, coli);
        rowinfile[coli] =  Mt(rowi, coli) + '0'; // forming string row before writing to file
     } 
 for(long ii=0; ii< Mt.cols() ; ii++){
     fileOUT << rowinfile[ii];
  }

     fileOUT << "\n";
  }
 
 fileIN.close();
 fileOUT.close();

} else {
   //  Situation 2 
   //  Block approach needed due to lack of memory

   if (!quiet ){
        message( " A block transpose is being performed due to lack of memory.  ");
        message( " Memory parameter availmemGb is set to " , max_memory_in_Gbytes , "gigabytes" );
        message( " If possible, increase availmemGb parameter. " );
    }

    // Calculate number of columns that can be read into available memory
    n_of_cols_to_be_read = max_mem_in_bytes * 1.0 / (3.5 * dims[0] * (bits_in_int/8.0)); //64 bit system
    // Calculate number of blocks needed
    long n_blocks = dims[1]/n_of_cols_to_be_read;
    if (dims[1] % n_of_cols_to_be_read != 0)
         n_blocks++;
    
    if (!quiet ) { 
         message( " Block Transpose of ASCII genotype file beginning ... " );
         message("  Due to marker data exceeding memory, data being processed in blocks. Number of blocks being processed is ", n_blocks);
    }
    // Block read and transpose - requires n_blocks passes through the 
    // ASCII input file which could be slow if file is large and memory low
    for(long b=0; b < n_blocks; b++){
         if (!quiet ) 
              message( " Processing block ... " , b , " of a total number of blocks of " , n_blocks );


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
     //    long counter = 0;
         if (!quiet ) {
              message("\n\n");
         }
         for(long rowi=0; rowi<dims[0]; rowi++){ 
               // read a line of data from ASCII file
              getline(fileIN, line);
              std::istringstream streamA(line);
  
             long coli=0 ;
             for(long ii=start_val; ii < end_val ; ii++){
                M(rowi, coli)  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
                coli++;
             }
        } // end for(rowi=0; rowi<dims[0]; rowi++)
       // transpose M


       MatrixXi Mt = M.transpose();


      // write out contents fo Mt to file (no spaces)a
 //     std::string rowinfile(Mt.cols(), '0');  // initialising string with 0's
// std::string* rowinfile = NULL;
// rowinfile = new string  [ Mt.cols() ];
std::vector<char> rowinfile( Mt.cols() ) ;


 for(long i=0; i < Mt.cols(); i++)
    rowinfile[i] = '0';


      for(long rowi=0; rowi<Mt.rows(); rowi++){
         for(long coli=0; coli<Mt.cols(); coli++){
             // fileOUT << Mt(rowi, coli);
            rowinfile[coli] =  Mt(rowi, coli) + '0'; // forming string row before writing to file
         }
         for(long ii=0; ii < Mt.cols(); ii++)
            fileOUT << rowinfile[ii]; // writing entire row of data
         fileOUT << "\n";
      }




      // return to the beginning of the input file
      fileIN.clear();
      fileIN.seekg(0, std::ios::beg);

     }     // end for blocks

fileIN.close();
fileOUT.close();


}  // end if else situation 


message(" The marker file has been Uploaded");


}  // end function 







//--------------------------------------------
// Calculation of transformed blup a values
//--------------------------------------------
// [[Rcpp::export]]
Eigen::MatrixXd calculate_reduced_a_rcpp ( CharacterVector f_name_ascii, double varG, 
                                           Eigen::Map<Eigen::MatrixXd> P,
                                           Eigen::Map<Eigen::MatrixXd>  y,
                                           double max_memory_in_Gbytes,  
                                           std::vector <long> dims,
                                           Rcpp::NumericVector  selected_loci,
                                           bool quiet,
                                           Rcpp::Function message)
{
  // function to calculate the BLUPs for the dimension reduced model. 
  // It is being performed in Rcpp because it makes use of Mt. 
  // Args
  // f_name_ascii    path + file name of Mt.bin
  // varG          variance of polygenic component
  // P             calculate in R
  // y             response/trait  but read in as a row matrix
  // max_memory_in_Gbytes  working memory in gigabytes
  // dims          dimension (row, column), of M.

std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

Eigen::MatrixXd
      ar(dims[1],1);  // column vector

std::ostringstream
      os;

Eigen::MatrixXd 
      nullmat = Eigen::MatrixXd::Zero(1,1);  // 0 matrix of size 1,1 for null purposes 


// const size_t bits_in_double = std::numeric_limits<double>::digits;


   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt))
double mem_bytes_needed =   ( dims[0]*dims[1] + dims[0]*dims[0] + dims[0] ) *  ( sizeof(double)/( 1000000000));

if (!quiet){
    // Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
    message("Inside internal function calculate_reduced_a_rcpp. Memory needed (gigabytes): ", mem_bytes_needed);
    message("Inside internal function calculate_reduced_a_rcpp. Memory available (gigabytes): ", max_memory_in_Gbytes);


}

if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory

   Eigen::MatrixXd
                   Mt;

   Mt = ReadBlock(fnamebin, 0, dims[0], dims[1]);
 //  Mt = ReadBlockFast(fnamebin, 0, dims[0], dims[1]);

  if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(long ii=0; ii < selected_loci.size() ; ii++)
          Mt.row(selected_loci(ii) ).setZero();
   }

   // ar  =    varG * Mt *  P   * y ;
    ar  =     P   * y ;
    ar  =    Mt * ar;
    ar  =    varG * ar;
} else {

      // calculation being processed in block form
      message(" Note:  Increasing availmemGb would improve performance... ");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algrebra but it is comes to the following
      long num_rows_in_block = (max_memory_in_Gbytes * 1000000000.0/sizeof(double) - dims[0] * dims[0] - dims[0])/dims[0] ;

    if (num_rows_in_block < 0){
        message("\n");
        message("Error:  availmemGb is set to " , max_memory_in_Gbytes );
        message("        Cannot even read in a single row of data into memory." );
        message("        Please increase availmemGb for this data set." );
        message("\n");
        message( "AM has terminated with errors\n" );
        return nullmat;

      }


      // blockwise multiplication

      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;

      if (!quiet ){
      message(" Block multiplication necessary. \n");
      message(" Number of blocks needing block multiplication is ... % d \n", num_blocks);
      } 
      for(long i=0; i < num_blocks; i++){
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[1])
            num_rows_in_block1 = dims[1] - start_row1;

          Eigen::MatrixXd
                  Mt;
            Mt = ReadBlock(fnamebin, start_row1, dims[0], num_rows_in_block1) ;
        //   Mt = ReadBlockFast(fnamebin, start_row1, dims[0], num_rows_in_block1) ;

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
                if(selected_loci(ii) >= start_row1 && selected_loci(ii) < start_row1 + num_rows_in_block1 )
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

       if (!quiet  )  message( "block done ... ");
      } // end for long




}  // end if mem_bytes_needed

  return(ar);

} // end function 




// internal function to remove a row from a dynamic matrix
void removeRow(Eigen::MatrixXd& matrix, unsigned long rowToRemove)
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
                                    Eigen::Map<Eigen::MatrixXd> inv_MMt_sqrt,
                                    Eigen::Map<Eigen::MatrixXd> dim_reduced_vara,
                                    double  max_memory_in_Gbytes,  
                                    std::vector <long> dims,
                                    Eigen::VectorXd  a,
                                    bool  quiet, 
                                    Rcpp::Function message)
{
// Purpose: to calculate the untransformed BLUP (a) values from the 
//          dimension reduced BLUP value estimates. 
//          It is necessary to have a block multiplication form of this function. 
//          Also, since the matrix multiplications are reliant upon the BLAS library, only 
//          double precision matrix multiplication is possible. This means, the Mt matrix must 
//          be converted into a double precision matrix which has a large memory cost.  
// Note:
//      1. dims is the row, column dimension of the Mt matrix


std::ostringstream
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

// const size_t bits_in_double = std::numeric_limits<double>::digits;
// const size_t bits_in_integer = std::numeric_limits<int>::digits;



   // Calculate memory footprint for Mt %*% inv(sqrt(MMt)) %*% var(a) %*% inv(sqrt(MMt)%*%M)
//AWG  double mem_bytes_needed =   ( dims[0]   +  2*dims[1]   + 1 ) *  (dims[1] * sizeof(double) /( 1000000000));
 double mem_bytes_needed =   ( 4   *dims[1]  *  dims[0] * sizeof(double))/1000000000;

if (!quiet){
//   Rprintf("Total memory (Gbytes) needed for a calculation is: %f \n",  mem_bytes_needed);
   message("Inside internal function calculate_a_and_vara_rcpp: Need memory (gigabytes)  ", mem_bytes_needed);
}





if(mem_bytes_needed < max_memory_in_Gbytes){
 // calculation will fit into memory
     Eigen::MatrixXd Mt = ReadBlock(fnamebin, 0, dims[1], dims[0]);
 //    Eigen::MatrixXd Mt = ReadBlockFast(fnamebin, 0, dims[1], dims[0]);


   if(!R_IsNA(selected_loci(0))){
   // setting columns to 0
   for(long ii=0; ii < selected_loci.size() ; ii++){
           Mt.row(selected_loci(ii)).setZero();
    }
   }




//   AWGans =    Mtd *  inv_MMt_sqrt  * a ;
//   AWGans =    Mt.cast<double>() *  inv_MMt_sqrt  * a ;
// std::clock_t    start;

//   start = std::clock();
    Eigen::MatrixXd  ans_part1 = inv_MMt_sqrt * a;


    ans.noalias() =   Mt  * ans_part1; 

   




// AWG    ans_part1.resize(0,0);  // erase matrix
   //  ans =    Mt *  inv_MMt_sqrt  * a ;
 //  Rprintf(" finished untransfomred BLUP values \n");



  // calculate untransformed variances of BLUP values
//  Eigen::MatrixXd var_ans_tmp_part1 =  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;

    Eigen::MatrixXd var_ans_tmp_part1 =   dim_reduced_vara * inv_MMt_sqrt;

    var_ans_tmp_part1 = inv_MMt_sqrt * var_ans_tmp_part1;



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
      message(" Increasing maxmemGb would improve performance... \n");

      // calculate the maximum number of rows in Mt that can be contained in the
      // block multiplication. This involves a bit of algebra but it is comes to the following
      // Strictly, 2 should be used here but I want some extra memory to play with 
      long num_rows_in_block =  max_memory_in_Gbytes * ( 1000000000) /
                             ( 4  *dims[1] *  sizeof(double) ) ;


      if (num_rows_in_block < 0){
        message("\n");
        message( "Error:  availmemGb is set to " , max_memory_in_Gbytes );
        message( "        Cannot even read in a single row of data into memory." );
        message( "        Please increase availmemGb for this data set." );
        message("\n");
        message(" multiple_locus_am has terminated with errors\n" );
        
       return Rcpp::List::create(Rcpp::Named("a")=0,
                            Rcpp::Named("vara") = 0);

      }
     


      // blockwise multiplication

      // find out the number of blocks needed
      long num_blocks = dims[0]/num_rows_in_block;
      if (dims[0] % num_rows_in_block)
                 num_blocks++;
      if (!quiet  ){
      message(" Block multiplication necessary. \n");
      message(" Number of blocks needing block multiplication is ... % d \n", num_blocks);
      } 
      for(long i=0; i < num_blocks; i++){
         message("Performing block iteration ... " , i );
         long start_row1 = i * num_rows_in_block;
         long num_rows_in_block1 = num_rows_in_block;
         if ((start_row1 + num_rows_in_block1) > dims[0])
            num_rows_in_block1 = dims[0] - start_row1;


           Eigen::MatrixXd Mt = ReadBlock(fnamebin, start_row1, dims[1], num_rows_in_block1) ;
       //   Eigen::MatrixXd Mt = ReadBlockFast(fnamebin, start_row1, dims[1], num_rows_in_block1) ;



        Eigen::MatrixXd
              vt1,
              ans_tmp1;
          
        Eigen::MatrixXd   var_ans_tmp(num_rows_in_block1,1);

            if(!R_IsNA(selected_loci(0))){
            // setting columns (or row when Mt) to 0
               for(long ii=0; ii < selected_loci.size() ; ii++)
               {
               // since we are now dealing with Mt, and blocking on columns, 
               // because columns are rows in Mt, then we have to be careful
               // that we do not select loci outside the block bounds. Also 
               // the values have to be adjusted based on the block number
                   if(selected_loci(ii) >= start_row1 && selected_loci(ii) < start_row1 + num_rows_in_block1 )
                   {   // selected loci index is in block 
                   long block_selected_loci = selected_loci(ii) - start_row1;
                   Mt.row(block_selected_loci).setZero();
                   }
                }   
            }
           //  ans_tmp  =  Mtd *  inv_MMt_sqrt  * a ;
             ans_tmp.noalias()  =   inv_MMt_sqrt  * a ;
             ans_tmp = Mt * ans_tmp;

            // variance calculation
            // vt.noalias() =  Mtd *  inv_MMt_sqrt * dim_reduced_vara * inv_MMt_sqrt;
             vt1.noalias() =  dim_reduced_vara * inv_MMt_sqrt;
             vt1           =  inv_MMt_sqrt * vt1;



        // performing quadratic form, remembering only diag elements are needed for variances. 
//        Rcout << " Performing ... Mt.row(j) * vt1 * ((Mt.row(j)).transpose())  " << endl;
//        for(long j=0; j < num_rows_in_block1; j++){
//           var_ans_tmp(j,0) =  Mt.row(j) * vt1 * ((Mt.row(j)).transpose()) ;
//        }
//        Rcout << "end of computing variances ... " << endl;
            // vt.noalias() =  Mt *  vt1;
           Eigen::MatrixXd vt;
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
    
             if (!quiet  )  message( "block done ... " );


      } // end for long



}  //  end if block update


  return Rcpp::List::create(Rcpp::Named("a")=ans,
                            Rcpp::Named("vara") = var_ans);


}






// [[Rcpp::export]]
bool  createM_ASCII_rcpp(CharacterVector f_name, CharacterVector f_name_ascii, 
                  CharacterVector  type,
                  std::string AA,
                  std::string AB, 
                  std::string BB,
                  double  max_memory_in_Gbytes,  std::vector <long> dims,
                  bool quiet,
                  Rcpp::Function message, 
                  std::string missing) 
{
  // Rcpp function to create space-removed ASCII file from ASCII and PLINK input files




// size_t found;


std::string 
   line; 


std::ofstream
   fileOUT;

//int 
//   genoval;

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
       CreateASCIInospace_PLINK(fname, fnameascii, dims, quiet, message);

   }  else {
      //-------------------------------------------
      // convert text file into ASCII file
      //-----------------------------------------
      // Here, we do not need to worry about the amount of memory because 
      // we are processing a line of the file at a time. This is not the case when 
      // creating a ASCII Mt because we have to read in blocks before we can 
      // transpose. 
      if (!quiet  )
          message(" A text file is being assumed as the input data file type. ");

       if ( 1.0 * memory_needed_in_Gb   > max_memory_in_Gbytes){
           bool it_worked = CreateASCIInospace(fname, fnameascii, dims, AA, AB, BB, quiet, message, missing);
           if (!it_worked) // an error has occurred in forming ascii file
                 return false;
       } else {
        //    bool it_worked =  CreateASCIInospaceFast(fname, fnameascii, dims, AA, AB, BB, quiet, message, missing);
            bool it_worked =  CreateASCIInospace(fname, fnameascii, dims, AA, AB, BB, quiet, message, missing);
          if (!it_worked) // an error has occurred in forming ascii file
                 return false;
       } 


   }  // end if type == "PLINK" 

//--------------------------------------
// Summary of Genotype File
//--------------------------------------

message( "\n\n                    Summary of Marker File  " );
message( "                   ~~~~~~~~~~~~~~~~~~~~~~~~   " );
message( " File type:                " , type  );
message(" File name:                " , fname );
message(" New ASCII file name:  " , fnameascii  );
message(" Number of individuals:    "     , dims[0] );
if (ftype == "PLINK"  ){
message(" Number of loci:           "  , (dims[1] -6)/2.0   );
} else {
message(" Number of loci:           "  , dims[1] );
}
message( " File size (gigabytes):       "  , memory_needed_in_Gb );
message(" Available memory (gigabytes):" , max_memory_in_Gbytes  );
message("\n\n" );


  return true; 


}











// [[Rcpp::export]]
Eigen::VectorXi  extract_geno_rcpp(CharacterVector f_name_ascii, 
                                   double  max_memory_in_Gbytes, 
                                    long  selected_locus, 
                                    std::vector<long> dims)
{
  std::string
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

  long 
     nind;

  nind = dims[0];




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
  //  Eigen::MatrixXd genoMat =  ReadBlockFast(fnamebin,  0, dims[1], dims[0]);

   column_of_genos = genoMat.col(selected_locus).cast<int>() ;
   


}  else {
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
              // Eigen::MatrixXd    
             //    genoMat_block1 ( ReadBlockFast(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;

              // removing rows that correspond to individuals with no 
              // trait data
           //   long start = start_row1;
           //   long  finish = start_row1 + num_rows_in_block1;


              // dealing with assigning column_of_genos when some values 
              // may be missing due to having been removed. 
              long colindx = start_row1;
              for(long j=start_row1; j< start_row1+num_rows_in_block1 ; j++){
                  column_of_genos(colindx) = genoMat_block1.col(selected_locus)(j-start_row1);
                  colindx++;

              } // end for j

          } // end for  i


} // end if max_memory

return(column_of_genos);

}






// [[Rcpp::export]]
Eigen::MatrixXd  calculateMMt_rcpp(CharacterVector f_name_ascii, 
                                   double  max_memory_in_Gbytes, int num_cores,
                                   Rcpp::NumericVector  selected_loci , std::vector<long> dims, 
                                   bool  quiet, Rcpp::Function message) 
{
// set multiple cores
Eigen::initParallel();
omp_set_num_threads(num_cores);
Eigen::setNbThreads(num_cores);
message(" Number of cores being used for calculation is .. ", num_cores);


std::string 
   line; 


std::ofstream
   fileOUT;

//int 
//   genoval;

std::string 
     fnamebin = Rcpp::as<std::string>(f_name_ascii);

// gpu will only work with double precision matrices in Eigen. 
// Had to change code to be double precision. 
Eigen::MatrixXd
    MMt(Eigen::MatrixXd(dims[0], dims[0]).setZero());





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
  //   Eigen::MatrixXd genoMat = ReadBlockFast(fnamebin,  0, dims[1], dims[0] );
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
    double part1 = -2 *  dims[1];
    double part2 = 4*dims[1] * dims[1] +  4 * max_memory_in_Gbytes  * 1000000000.0/sizeof(double);
    part2 = sqrt(part2);
    long num_rows_in_block = (part1 + part2)/2.2;
    message( "number of rows in block is " , num_rows_in_block);

           // blockwise multiplication

          // find out the number of blocks needed
        message( " 1 ");
          long num_blocks = dims[0]/num_rows_in_block;



          if (dims[0] % num_rows_in_block)
                 num_blocks++;


          for(long i=0; i < num_blocks; i++){
              long start_row1 = i * num_rows_in_block;
              long num_rows_in_block1 = num_rows_in_block;
              if ((start_row1 + num_rows_in_block1) > dims[0])
                     num_rows_in_block1 = dims[0] - start_row1;


                Eigen::MatrixXd    
                     genoMat_block1 ( ReadBlock(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;
        //       Eigen::MatrixXd    
        //            genoMat_block1 ( ReadBlockFast(fnamebin,  start_row1, dims[1], num_rows_in_block1)) ;
              Eigen::MatrixXd    
                   MMtsub(Eigen::MatrixXd(num_rows_in_block1, num_rows_in_block1).setZero());

             if(!R_IsNA(selected_loci(0) )){
             // setting columns to 0
             for(long ii=0; ii < selected_loci.size() ; ii++)
                genoMat_block1.col(selected_loci(ii)).setZero();
             }
             // Rcpp::Rcout << "  Block 1  "  << std::endl;
             // Rcpp::Rcout << genoMat_block1.rows() << std::endl;
             // Rcpp::Rcout << genoMat_block1.cols() << std::endl;
              MMtsub.noalias() = genoMat_block1 * genoMat_block1.transpose(); 
              //          i            j            num rows               num   cols
              MMt.block(start_row1, start_row1, num_rows_in_block1, num_rows_in_block1) = MMtsub;
       
              for(long j=i+1;j<num_blocks; j++){
                   long start_row2 = j * num_rows_in_block;
                   long num_rows_in_block2 = num_rows_in_block;
                   if ((start_row2 + num_rows_in_block2) > dims[0])
                          num_rows_in_block2 = dims[0] - start_row2;
                     Eigen::MatrixXd    
                        genoMat_block2 ( ReadBlock(fnamebin,  start_row2, dims[1], num_rows_in_block2)) ;

              //       Eigen::MatrixXd    
              //          genoMat_block2 ( ReadBlockFast(fnamebin,  start_row2, dims[1], num_rows_in_block2)) ;



                   Eigen::MatrixXd    MMtsub(Eigen::MatrixXd(num_rows_in_block1, num_rows_in_block2).setZero());

                  if(!R_IsNA(selected_loci(0) )){
                   // setting columns to 0
                   for(long jj=0; jj < selected_loci.size() ; jj++)
                      genoMat_block2.col(selected_loci(jj)).setZero();
                   }
            //  Rcpp::Rcout << " Block 2 " << std::endl;
            //  Rcpp::Rcout << genoMat_block2.rows() << std::endl;
            //  Rcpp::Rcout << genoMat_block2.cols() << std::endl;
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





unsigned long getNumColumns(std::string fname) {
	unsigned long numcols = 0;

	std::string line;

	std::ifstream fileIN(fname.c_str());

	getline(fileIN, line);
	std::istringstream streamA(line);

	std::string token;

	numcols = line.length();

	return numcols;

}

unsigned long getNumRows(std::string fname) {
	unsigned long numrows = 0;

	std::string line;

	std::ifstream fileIN(fname.c_str());

	std::istringstream streamA(line);

	// Determine number of columns in file
	fileIN.clear(); // returns to beginning of line
	fileIN.seekg(0, std::ios::beg);

	// Determine number of rows in file
		while(fileIN.good()){
			while(getline(fileIN, line )){
				numrows++;
			}
		}

	return numrows;

}


