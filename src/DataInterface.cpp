//
// Created by Mohamed Sneiba on 9/28/18.
//
/**
 * @file DataInterface.cpp
 *
 * Implementation of the DataInterface class
 *
 */
#include "DataInterface.h"
#include <armadillo>
#include <string>
#include <iostream>
#include <vector>

using namespace arma;
using namespace std;

/**
* Main Constructor of the DataInterface Class, it assigns default values to the three object properties
* @param input contains the path to the input file. By default it contains the value '/resources/in.csv'.
* @param output contains the path to the output file where the results will be store. By default it contains the value '/resources/out.csv'.
* @param delim contains the delimiter of the data inside the input and the output file. By default the delimiter is ','.
*/
DataInterface::DataInterface(string input, string output, string delim) :
	inputFileName(input), outputFileName(output), delimiter(delim)
{

}

/**
* This function allows us to split by a delimiter a given string into a vector of string.
* It will be useful for us in order the extract the data from each line in the input file.
* @param stringToBeSplitted represents the string that contains all the data and need to be splitted.
* @param delimiter represents the delimiter to use when splitting the given string.
* @return This function returns a vector of string contains all the data that was separated by the delimiter before.
*/
vector<string> split(string stringToBeSplitted, string delimiter)
{
	vector<string> splittedString;
	int startIndex = 0;
	int  endIndex = 0;
	while((endIndex = stringToBeSplitted.find(delimiter, startIndex)) < stringToBeSplitted.size())
	{

		string val = stringToBeSplitted.substr(startIndex, endIndex - startIndex);
		splittedString.push_back(val);
		startIndex = endIndex + delimiter.size();

	}
	if(startIndex < stringToBeSplitted.size())
	{
		string val = stringToBeSplitted.substr(startIndex);
		splittedString.push_back(val);
	}
	return splittedString;
}

/** Getter of the inputFileName property
 *
 * @return This function returns the inputFileName value
 *
 */
string DataInterface::getInputFileName()
{
	return this->inputFileName;
}
/** Getter of the outputFileName property
 *
 * @return This function returns the outputFileName value
 *
 */
string DataInterface::getOutputFileName()
{
	return this->outputFileName;
}
/** Getter of the delimiter property
 *
 * @return This function returns the delimiter value
 *
 */
string DataInterface::getDelimiter()
{
	return this->delimiter;
}

/** Setter of the inputFileName property
 *
 * This function sets the value of inputFileName
 *
 * @param input
 *
 */
void DataInterface::setInputFileName(string input)
{
	this->inputFileName = input;
}
/** Setter of the outputFileName property
 *
 * This function sets the value of outputFileName
 *
 * @param output
 *
 */
void DataInterface::setOutputFileName(string output)
{
	this->outputFileName = output;
}
/** Setter of the delimiter property
 *
 * This function sets the value of delimiter
 *
 * @param delim
 *
 */
void DataInterface::setDelimiter(string delim)
{
	this->delimiter = delim;
}

/**
* This method allows us to read the data from the input file
* @param n This parameter is passed by address so we can put in it the value of n_max that we just read from the input file
* @return This method returns a vector containing the values of z that we just read from the input file.
*/

rowvec DataInterface::read(int& n)
{
	ifstream in(this->inputFileName);
	if(in)
	{
		string line = "";
		vector<string> data;
		vector<double> z;
		getline(in, line);
		getline(in, line);
		data = split(line, this->delimiter);
		n = stoi(data[0]);

		for(int i = 1; i < data.size(); i++)
		{
			z.push_back(stod(data[i]));
		}

		return rowvec(z);

	}
	else
	{
		cout << "ERROR While opening the file." << endl;
		return NULL;
	}
}

/**
* This method allows us to write the result of the computations we have made into the output file. This will allow us later on to plot the data.
* @param z This parameter represents the inital z vector, the function will write it in the first column of the output file
* @param w This parameter is the matrix of all the \f$\psi\f$ that we computed and write it in the ouput file.
*/
void DataInterface::write(rowvec z, mat w)
{
	ofstream out(this->outputFileName);
	if(out)
	{
		string line = "z";
		mat wt = w.t();
		for(int i = 0; i < wt.n_cols; i++)
		{
			line = line + ",w[" + to_string(i) + "]";
		}
		out << line << endl;
		for(int i = 0; i < wt.n_rows; i++)
		{
			line = to_string(z[i]) + "," + to_string(wt.at(i, 0));
			for(int j = 1; j < wt.n_cols; j++)
			{
				line = line + "," + to_string(wt.at(i, j));
			}

			out << line << endl;
		}

	}
	else
	{
		cout << "ERROR While opening the file." << endl;

	}

}
/// Destructor of the Object DataInterface
DataInterface::~DataInterface()
{

}
