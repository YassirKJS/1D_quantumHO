/**
 * @file DataInterface.h
 *
 * Header of the DataInterface class
 *
 */
#ifndef DEF_DATA_INTERFACE
#define DEF_DATA_INTERFACE

#include <armadillo>
#include <string>

using namespace arma;
using namespace std;

/**
 * @class DataInterface
 *
 * This class helps us get the entry data from a .csv file and write the result of the computations in a .csv output file
 *
 */

class DataInterface
{
	private:
		/// the property that stores the path of the file that contains the entry data
		string inputFileName;
		/// the property that stores the path of the file where to write the results.
		string outputFileName;
		/// the data delimiter in each line.
		string delimiter;
	public :

		/// Main Constructor of the DataInterface Class, it assigns default values to the three object properties
		DataInterface(string input = "../resources/in.csv", string output = "../resources/out.csv", string delim = ",");

		/// This method reads the content of the input file and returs the vector z and the n_max that is passed by adress in the function's parameters.
		rowvec read(int&);
		/// This method writes the content of the vector z, and the matrix containing all the values of \f[\psi\f]
		void write(rowvec, mat);

		/// The getter of the inputFileName object property
		string getInputFileName();
		/// The getter of the outputFileName object property
		string getOutputFileName();
		/// The getter of the delimiter object property
		string getDelimiter();

		/// The setter of the inputFileName object property
		void setInputFileName(string input);
		/// The setter of the outputFileName object property
		void setOutputFileName(string output);
		/// The setter of the delimiter object property
		void setDelimiter(string delim);

		/// Destructor of the Object DataInterface
		~DataInterface();
};

#endif