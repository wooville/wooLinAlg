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
void printMatrix(wooMatrix2D<T> matrix)
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
void printVector(wooVector<T> inputVector)
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
		
		printMatrix(testMatrix);
		
		cout << endl;
		cout << "Computing eigenvector and eigenvalue..." << endl;
		double eigenValue;
		wooVector<double> eigenVector;
		wooEIG_PIt<double>(testMatrix, eigenValue, eigenVector);
		
		cout << "Eigenvector: " << endl;
		printVector(eigenVector);
		cout << "Eigenvalue = " << eigenValue << "." << endl;
		cout << endl;
	}

    cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "Inverse-Power Iteration Method." << endl;
	cout << "**********************************************" << endl;
	cout << endl;	
	
	{
		cout << "Testing with a simple 3x3 matrix:" << endl;
		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		printMatrix(testMatrix);
		cout << endl;
		
		// Define the eigenvalues.
		std::vector<double> eigenValues = {1.5962551, -0.37253087, 0.02627577};
		
		// Setup a vector for the eigenvector.
		wooVector<double> eigenVector(3);
		
		// Loop through each eigenvalue.
		for (auto currentValue : eigenValues)
		{
			cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
			int returnStatus = wooEIG_InvPIt<double>(testMatrix, currentValue, eigenVector);
			printVector(eigenVector);
			
			if (returnStatus == WOOEIG_MAXITERATIONSEXCEEDED)
				cout << "*** Maximum iterations exceeded ***" << endl;
			
			cout << endl;
		}
		
	}
	
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "New isSymmetric() function in wooMatrix2D class." << endl;
	cout << "**********************************************" << endl;
	cout << endl;			
	
	{
		cout << "Testing with a simple symmetric matrix." << endl;
		std::vector<double> simpleData = {4.0, -7.0, 6.0, -7.0, -2.0, 13.0, 6.0, 13.0, 5.0};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		printMatrix(testMatrix);
		cout << endl;
		
		cout << "Matrix is symmetric: ";
		if (testMatrix.isSymmetric())
			cout << "True." << endl;
		else
			cout << "False." << endl;
	}
	
	{
		cout << "Testing with a simple non-symmetric matrix." << endl;
		std::vector<double> simpleData = {4.0, -7.0, 6.0, 7.0, -2.0, 13.0, 6.0, 13.0, 5.0};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		printMatrix(testMatrix);
		cout << endl;
		
		cout << "Matrix is symmetric: ";
		if (testMatrix.isSymmetric())
			cout << "True." << endl;
		else
			cout << "False." << endl;
	}		
	
	cout << endl << endl;	
	cout << "**********************************************" << endl;
	cout << "Testing eigenvalue and eigenvector code." << endl;
	cout << "Eigenvalues by QR decomposition." << endl;
	cout << "**********************************************" << endl;
	cout << endl;		
	
	{
		cout << "Testing with a simple (non-symmetric) 3x3 matrix:" << endl;
		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		printMatrix(testMatrix);
		cout << endl;
		
		// Compute the eigenvalues.
		std::vector<double> eigenValues;
		int returnStatus = wooEIG_QR(testMatrix, eigenValues);
		
		if (returnStatus == WOOEIG_MAXITERATIONSEXCEEDED)
			cout << ">>> Maximum iterations exceeded <<<" << endl;
			
		if (returnStatus == WOOEIG_MATRIXNOTSYMMETRIC)
			cout << ">>> Matrix not symmetric. <<<" << endl;	
		
		// Display the eigenvalues.
		cout << "The estimated eigenvalues are:" << endl;
		for (auto currentValue : eigenValues)
			cout << std::setprecision(6) << currentValue << " ";
			
		cout << endl << endl;	
	}
	
	{
		cout << "Testing with a simple (symmetric) 3x3 matrix:" << endl;
		std::vector<double> simpleData = {6.0, 5.5, -1.0, 5.5, 1.0, -2.0, -1.0, -2.0, -3.0};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		printMatrix(testMatrix);
		cout << endl;
		
		// Compute the eigenvalues.
		std::vector<double> eigenValues;
		int returnStatus = wooEIG_QR(testMatrix, eigenValues);
		
		if (returnStatus == WOOEIG_MAXITERATIONSEXCEEDED)
			cout << ">>> Maximum iterations exceeded <<<" << endl;
			
		if (returnStatus == WOOEIG_MATRIXNOTSYMMETRIC)
			cout << ">>> Matrix not symmetric. <<<" << endl;	
		
		// Display the eigenvalues.
		cout << "The estimated eigenvalues are:" << endl;
		for (auto currentValue : eigenValues)
			cout << std::setprecision(6) << currentValue << " ";
			
		cout << endl << endl;
			
		// Setup a vector for the eigenvector.
		wooVector<double> eigenVector(3);
		
		// Loop through each eigenvalue.
		for (auto currentValue : eigenValues)
		{
			cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
			int returnStatus = wooEIG_InvPIt<double>(testMatrix, currentValue, eigenVector);
			printVector(eigenVector);
			
			if (returnStatus == WOOEIG_MAXITERATIONSEXCEEDED)
				cout << "*** Maximum iterations exceeded ***" << endl;
			
			cout << endl;
		}			
			
		cout << endl << endl;	
	}	
	
	{
		cout << "Testing with an example that should have complex eigenvalues:" << endl;
		std::vector<double> simpleData = {4.0, -6.0, 8.0, 7.0, 9.0, -5.0, 9.0, -6.0, -4.0};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		printMatrix(testMatrix);
		cout << endl;
		
		// Compute the eigenvalues.
		std::vector<double> eigenValues;
		int returnStatus = wooEIG_QR(testMatrix, eigenValues);
		
		if (returnStatus == WOOEIG_MATRIXNOTSYMMETRIC)
			cout << ">>> Matrix not symmetric. <<<" << endl;		
		
		if (returnStatus == WOOEIG_MAXITERATIONSEXCEEDED)
			cout << ">>> Maximum iterations exceeded <<<" << endl;			
		
		// Display the eigenvalues.
		cout << "The estimated eigenvalues are:" << endl;
		for (auto currentValue : eigenValues)
			cout << std::setprecision(6) << currentValue << " ";
			
		cout << endl << endl;	
		
		// Setup a vector for the eigenvector.
		wooVector<double> eigenVector(3);		
		
		// Loop through each eigenvalue.
		for (auto currentValue : eigenValues)
		{
			cout << "Estimated eigenvector for eigenvalue " << currentValue << " = " << endl;
			int returnStatus = wooEIG_InvPIt<double>(testMatrix, currentValue, eigenVector);
			printVector(eigenVector);
			
			if (returnStatus == WOOEIG_MAXITERATIONSEXCEEDED)
				cout << "*** Maximum iterations exceeded ***" << endl;
			
			cout << endl;
		}
		
		cout << endl << endl;
	}
}