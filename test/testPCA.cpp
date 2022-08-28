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
#include "../wooPCA.h"

using namespace std;

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing principal component analysis code." << endl;
	cout << "**********************************************" << endl;
	cout << endl;
	
	{
		cout << "Testing with 100 observations of 3 variables:" << endl;
		
		// Read data from the .CSV file.
		string rowData;
		string number;
		stringstream rowDataStream;
		std::vector<double> testData;
		int numRows = 0;
		int numCols = 0;
		ifstream inputFile("wooPCATestData.csv");
		
		/* If the file open successfully then do stuff. */
		if (inputFile.is_open())
		{		
			cout << "Opened file successfully..." << endl;
			
			while (!inputFile.eof())
			{
			
				// Read the next line.
				getline(inputFile, rowData);
				
				// Loop through and extract the individual numbers.
				rowDataStream.clear();
				rowDataStream.str(rowData);
				
				if (numRows < 1)
					numCols = 0;
				
				while (rowDataStream.good())
				{
					getline(rowDataStream, number, ',');
					testData.push_back(atof(number.c_str()));
					
					if (numRows < 1)
						numCols++;
				}
				numRows++;
			}
						
			// Close the file.
			inputFile.close();
			
			//numRows--;
			//testData.pop_back();
            
			cout << "completed reading file..." << endl;
			cout << "Read " << numRows << " observations of " << numCols << " variables." << endl;
			cout << "Constituting " << testData.size() << " elements in total." << endl;
			
			// Form into a matrix.
			wooMatrix2D<double> X (numRows, numCols, testData);
			
			// compute the covariance matrix.
			std::vector<double> columnMeans = wooPCA::computeColMeans(X);
			wooMatrix2D<double> X2 = X;
			wooPCA::subtractColumnMeans(X2, columnMeans);
			
			wooMatrix2D<double> covX = wooPCA::computeCovariance(X2);
			cout << endl;
			cout << "Giving the covariance matrix as: " << endl;
			covX.printMatrix();
			
			// compute the eigenvectors.
			wooMatrix2D<double> eigenvectors;
			int testResult = wooPCA::computeEigenvectors(covX, eigenvectors);
			cout << endl;
			cout << "And the eigenvectors as: " << endl;
			eigenvectors.printMatrix();
			
			// Test the overall function.
			cout << endl;
			cout << "Testing overall function..." << endl;
			wooMatrix2D<double> eigenvectors2;
			int testResult2 = wooPCA::wooPCA(X, eigenvectors2);
			cout << "testResult2 = " << testResult2 << endl;
			cout << "And the final eigenvectors are:" << endl;
			eigenvectors2.printMatrix();
			
			// Test dimensionality reduction.
			cout << endl;
			cout << "Testing dimensionality reduction." << endl;
			cout << "Starting with X which has " << X.getNumRows() << " rows and " << X.getNumCols() << " columns." << endl;
			cout << endl;
			cout << "Using only the first two principal components:" << endl;
			wooMatrix2D<double> V, part2;
			eigenvectors.separate(V, part2, 2);
			V.printMatrix(8);
			cout << endl;
			
			wooMatrix2D<double> newX = (V.transpose() * X.transpose()).transpose();
			cout << "Result has " << newX.getNumRows() << " rows and " << newX.getNumCols() << " columns." << endl;
			
			// Open a file for writing
			ofstream outputFile("wooPCATestData_Reduced.csv");
			if (outputFile.is_open())
			{
				for (int i=0; i<newX.getNumRows(); ++i)
				{
					outputFile << newX.getElement(i, 0) << "," << newX.getElement(i, 1) << endl;
				}
				outputFile.close();
			}
			
		}
	}
}