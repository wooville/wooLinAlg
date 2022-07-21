#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <random>

#include "../wooMatrix2D.h"
#include "../wooVector.h"
#include "../wooLinSolve.h"

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
  cout << std::fixed << std::setprecision(3) << inputVector.getElement(row) << endl;
	}    
}

int main()
{
	cout << "Code to test wooMatrix2D" << endl;
	cout << "Testing conversion of matrix to row echelon form." << endl;
	cout << endl;
	
	// Generate a matrix to test things with.
    std::vector<double> simpleData = {1.0, 3.0, -1.0, 13.0, 4.0, -1.0, 1.0, 9.0, 2.0, 4.0, 3.0, -6.0};
    wooMatrix2D<double> testMatrix(3, 4, &simpleData);
  
  cout << "Original matrix:" << endl;
  PrintMatrix(testMatrix);
  cout << endl;
  
  // Convert to row echelon form.
  wooMatrix2D<double> rowEchelonMatrix = testMatrix.rowEchelon();
  
  cout << "Converted to row echelon form:" << endl;
  PrintMatrix(rowEchelonMatrix);
  cout << endl;
  
  // Define another matrix as the first part of our system of linear equations.
  std::vector<double> simpleData2 = {1.0, 3.0, -1.0, 4.0, -1.0, 1.0, 2.0, 4.0, 3.0};
  wooMatrix2D<double> aMat(3, 3, &simpleData2);
  cout << "We setup the equations in the form of Ax = b, where A = " << endl;
  PrintMatrix(aMat);
  cout << endl;
  
  // Define a vector to hold the RHS of our system of linear equations.
  std::vector<double> vectorData {13.0, 9.0, -6.0};
  wooVector<double> bVec {vectorData};
  cout << "And b = " << endl;
  PrintVector(bVec);
  cout << endl;
  
  // Call the wooLinSolve function.
  wooVector<double> testResult = wooLinSolve<double>(&aMat, &bVec);
  cout << "And the final result is:" << endl;
  PrintVector(testResult);
  cout << endl;
  
  // ***************************************************************************************************
  // Try some random tests.
  std::random_device myRandomDevice;
  std::mt19937 myRandomGenerator(myRandomDevice());
  std::uniform_real_distribution<double> myDistribution(-25.0,25.0);
  
  int numUnknowns = 10;
  
  std::vector<double> coefficientData;
  std::vector<double> unknownData;
  // Populate the coefficient data.
  for (int i=0; i<(numUnknowns * numUnknowns); ++i)
  {
  	double randomNumber = myDistribution(myRandomGenerator);
  	coefficientData.push_back(randomNumber);
  }
  cout << "A random coefficient matrix = " << endl;
  wooMatrix2D<double> coefficientMatrix(numUnknowns, numUnknowns, &coefficientData);
  PrintMatrix(coefficientMatrix);
  cout << endl;
  
  cout << "And the random unknown values = " << endl;
  for (int i=0; i<numUnknowns; ++i)
  {
  	double randomNumber = myDistribution(myRandomGenerator);
  	unknownData.push_back(randomNumber);
  }
  wooVector<double> unknownVector {unknownData};
  PrintVector(unknownVector);
  cout << endl;
  
  cout << "Compute the equation results = " << endl;
  wooVector<double> systemResult = coefficientMatrix * unknownVector;
  PrintVector(systemResult);
  cout << endl;
  
  cout << "Attempt to solve the linear system..." << endl;
  wooVector<double> compSolution = wooLinSolve<double>(&coefficientMatrix, &systemResult);
  PrintVector(compSolution);
  cout << endl;
  
  cout << "And compare the actual result with the computed solution..." << endl;
  wooVector<double> errorVector = unknownVector - compSolution;
  PrintVector(errorVector);
  cout << endl;
  
	return 0;
}   