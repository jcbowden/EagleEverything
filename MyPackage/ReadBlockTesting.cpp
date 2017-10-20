#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <cstring>
#include <fcntl.h>
#include <time.h>
#include <bitset>

// Eigen package
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Eigen::MatrixXd  ReadBlock(std::string asciifname,
		long start_row,
		long numcols,
		long numrows_in_block);

Eigen::MatrixXd  ReadBlockFast(std::string asciifname,
		long start_row,
		long numcols,
		long numrows_in_block);

unsigned long getNumColumns(std::string fname);
unsigned long getNumRows(std::string fname);

char* mapFileFromDiscBlocked(const char * file_name, unsigned long &sizeUsed, unsigned long &sizeActual, unsigned long blockSize, unsigned long blockNum);

int main() {
	// File to be read
	string FILENAME = "/data/ste709/FilesForRyan/geno2.txt";

	// Number of rows to read per block iteration
	// used for testing purposes
	unsigned short inc = 3;
	bool quiet = false;

	// Get the number of columns in a row of data
	// Adapted from Andrew's code
	unsigned long numFileCols = getNumColumns(FILENAME.c_str());
	unsigned long numFileRows = getNumRows(FILENAME.c_str());


	for ( unsigned int i = 0; i < numFileRows - inc; i += inc)
	{
		// Read block of data using mmap() technique
		// Parameters: Filename, starting row, columns per line, number of lines to read
		MatrixXd tmp = ReadBlockFast(FILENAME, i, numFileCols, inc);
		MatrixXd tmp2 = ReadBlock(FILENAME, i, numFileCols, inc);

		// Use as checksum to check the results of the read process
		if ( tmp.sum() == tmp2.sum())
		{
			printf("Test passed for loading lines %d to %d\n", i, i+inc);
			if ( !quiet ) {
				printf("ReadBlockFast\tSum: %10f \n", tmp.sum());
				printf("ReadBlock\tSum: %10f \n", tmp2.sum());
			}
		}else {
			printf("Test failed for lines %d to %d\n", i, i+inc);
			exit(1);
		}
	}

	return 0;
}

// Ryan's ReadBlock code which uses mmap() system call
// Not fully tested, use with caution
Eigen::MatrixXd  ReadBlockFast(std::string asciifname,
		long start_row,
		long numcols,
		long numrows_in_block) {

	// Start settable parameters
	const unsigned long long maxMemory = 16ull*1024*1024*1024; // Mbytes
	// End settable parameters

	int pagesize = getpagesize();
	unsigned long sizeUsed = 0;
	unsigned long sizeActual = 0;

	// Used for offseting to a newline from a loaded block
	unsigned long realPos = (start_row)*(numcols+1);
	unsigned long allignedPos = realPos - (realPos % pagesize);
	unsigned long offsetCol = (realPos - allignedPos) % numcols;
	unsigned long offsetRow = floor((realPos - allignedPos) / numcols);

	// Debugging information
	/*cout << "NUM COLS: " << numcols << endl;
	/cout << "NUM ROWS: " << numrows_in_block << endl;
	cout << "REAL: " << realPos << endl;
	cout << "ALLIGNED: " << allignedPos << endl;
	cout << "Offset Col: " << offsetCol << endl;
	cout << "Offset Row: " << offsetRow << endl;
	cout << "Eigen Dimensions: (" << numrows_in_block << ", " << numcols << ")" << endl;
	cout << "Starting Point: " << offsetRow*(numcols) + offsetCol << endl;
	cout << "To Read: " << numrows_in_block * (numcols+1) << endl;
	*/

	// Eigen matrix to store block of data in
	Eigen::MatrixXd M(numrows_in_block, numcols) ;

	// Offset point to start reading mapped file from
	unsigned long startingPoint = offsetRow*(numcols) + offsetCol;

	// Memory required to read block of data
	unsigned long long requiredMemory = startingPoint + ( numrows_in_block*(numcols+1)) * 8; // bytes

	// Memory usage check
	if (requiredMemory <= maxMemory) {

		// Memory map file
		char* dataFile = mapFileFromDiscBlocked(asciifname.c_str(), sizeUsed, sizeActual, requiredMemory, allignedPos);

		// Used to load through the mapped file
		unsigned long rowInc = 0;
		unsigned long colInc = 0;

		// Loop through file data, reads data in as column major (Default for Eigen)
		for (unsigned int i = 0; i < numrows_in_block * (numcols+1); i++) {
			// If newline character is encountered
			if ( dataFile[startingPoint+i] == '\n') {
				colInc++;
				rowInc = 0;
			} else {
				// Read value from data file and convert to number format and subtract one
				// 0 -> -1, 1 -> 0, 2 -> 1 as per original ReadBlock code
				signed int value = (dataFile[startingPoint+i]  - '0') - 1;

				// Store value in Eigen matrix
				M(colInc,rowInc) = value;

				// Used for debugging purposes
				// cout << "Row: " << rowInc << " :: " << "Col: " << colInc << " :: " << value << endl;

				rowInc++;
			}

		}

		// Done using memory mapped file, release it
		munmap(dataFile, sizeUsed);

	} else {
		cout << "ERROR: Insufficent memory allocated to reading block of data." << endl;
		cout << "Required: " << requiredMemory/1024/1024/1024 << endl;
		exit(1);
	}

	return M;
}

/*
Eigen::Map<MatrixXd>  ReadBlock(std::string binfname,
		long start_row,
		long numcols,
		long numrows_in_block) {

	// STUB
	double result[10];
	Map<MatrixXd> tmp = Map<MatrixXd>(result, 5, 2);
	return tmp;
}
 */

char* mapFileFromDiscBlocked(const char * file_name, unsigned long &sizeUsed, unsigned long &sizeActual, unsigned long blockSize, unsigned long offset) {
	// If set to true debugging information is displayed
	// Otherwise these messages are suppressed
	bool debugMsgs = true;

	// Read file size and system page file size
	int pagesize = getpagesize();

	// Open file descriptor for given file name
	// Read-only permission
	int fd;

	// Open file and get file descriptor
	try {
		fd = open (file_name, O_RDONLY);
		// File descriptor returned as -1 if file cannot be opened
		if ( -1 == fd) throw "File could not be opened for reading";
	}
	catch (const char* msg) {
		cout << "ERROR: " << msg << endl;
		exit(-1);
	}

	// Round up file size to next multiple
	// of system page size;
	sizeActual = blockSize;
	sizeUsed = sizeActual + (pagesize - (sizeActual % pagesize));
	// if (debugMsgs) printf("Memory used to map file: %2.5f Mbytes\n", sizeUsed/1024/1024);

	// Memory mapped data file
	char *fileMemMap;

	// Map file to memory using mmap() system function
	// File may now be treated as a character array in memory
	// Syntax:
	// 		void *mmap(void *addr, size_t length,
	//			  int prot, int flags, int fd, off_t offset);
	// Permissions:
	//		PROT_READ: Pages in memory only allow read-only operations
	//		MAP_PRIVATE: Changes not visible to other processes
	//					 Underlying file is no altered

	// cout << "Page Size " << pagesize << endl;
	// cout << "Size Used " << sizeUsed << endl;
	// cout << "Size Actual " << sizeActual << endl;

	fileMemMap = (char *) mmap (0, sizeUsed, PROT_READ, MAP_PRIVATE, fd, offset);

	//if (madvise(fileMemMap, sizeUsed-1024, MADV_WILLNEED | MADV_SEQUENTIAL) == -1) {
	//	cout << "madvise error" << endl;
	//	return NULL;
	//}

	// Close the file descriptor
	// Frees resources associated with the file descriptor
	close(fd);

	return fileMemMap;
}

unsigned long getNumColumns(std::string fname) {
	unsigned long numcols = 0;

	string line;

	std::ifstream fileIN(fname.c_str());

	getline(fileIN, line);
	istringstream streamA(line);

	string token;

	numcols = line.length();

	return numcols;

}

unsigned long getNumRows(std::string fname) {
	unsigned long numrows = 0;

	string line;

	std::ifstream fileIN(fname.c_str());

	istringstream streamA(line);

	// Determine number of columns in file
	fileIN.clear(); // returns to beginning of line
	fileIN.seekg(0, ios::beg);

	// Determine number of rows in file
		while(fileIN.good()){
			while(getline(fileIN, line )){
				numrows++;
			}
		}

	return numrows;

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
		//Rcpp::stop(os.str() );
	}

	for( long rr=0; rr < (start_row + numrows_in_block) ; rr++){
		// read a line of data from ASCII file
		getline(fileIN, line);
		if(rr >= start_row){
			istringstream streamA(line);
			for( long ii=0; ii < numcols  ; ii++){
				int tmp  = line[ii] - '0'; // trick to removes ASCII character offset for numbers
				M(rowi, ii) = tmp - 1;   // converting data to -1, 0, 1
			}
			rowi++;
		} // end if rr
	} // end for(rr

	// Close the ascii file
	fileIN.close();


	return M;

}
