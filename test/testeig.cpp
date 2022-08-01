#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <fstream>

#include "../wooMatrix2D.h"
#include "../wooVector.h"
#include "../wooEIG.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void PrintMatrix(wooMatrix2D<T> matrix)
{
	int nRows = matrix.getNumRows();
	int nCols = matrix.getNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    cout << std::fixed << std::setprecision(3) << matrix.getElement(row, col) << "  ";
    }
	cout << endl;
	}    
}

// A simple function to print a vector to stdout.
template <class T>
void PrintVector(wooVector<T> inputVector)
{
	int nRows = inputVector.getNumDims();
	for (int row = 0; row<nRows; ++row)
  {
  cout << std::fixed << std::setprecision(6) << inputVector.getElement(row) << endl;
	}    
}

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "Power Iteration Method." << endl;
	cout << "**********************************************" << endl;
	cout << endl;
	
	{
		cout << "Testing with simple 3x3 matrix:" << endl;
		
		//std::vector<double> simpleData = {1.0, 2.0, 3.0, 4.0};
		//wooMatrix2D<double> testMatrix(2, 2, simpleData);
	
		std::vector<double> simpleData = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		
		PrintMatrix(testMatrix);
		
		cout << endl;
		cout << "Computing eigenvector and eigenvalue..." << endl;
		double eigenValue;
		wooVector<double> eigenVector;
		wooEIG_PIt<double>(testMatrix, eigenValue, eigenVector);
		
		cout << "Eigenvector: " << endl;
		PrintVector(eigenVector);
		cout << "Eigenvalue = " << eigenValue << "." << endl;
		cout << endl;
	}
}