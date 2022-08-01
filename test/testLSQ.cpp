/* *************************************************************************************************
	TestCode_wooLSQ
	
	  Code to test the linear least squares implementation in the wooLSQ function. Also tests matrix 
	 transpose.
	
	*** INPUTS ***
	
	None
															
	*** OUTPUTS ***
	
	INT				Flag indicating success or failure of the process.
	Created as part of the wooLinAlg linear algebra library, which is intended to be primarily for
	educational purposes. For more details, see the corresponding videos on the QuantitativeBytes
	YouTube channel at:
	
	www.youtube.com/c/QuantitativeBytes								
	************************************************************************************************* */

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
#include "../wooLSQ.h"

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
	cout << "Code to test linear least squares implementation" << endl;
	cout << "Testing matrix transpose." << endl;
	cout << endl;
	
	// Test matrix transpose.
	{
	  std::vector<double> simpleData = {1.0, 3.0, -1.0, 13.0, 4.0, -1.0, 1.0, 9.0, 2.0, 4.0, 3.0, -6.0};
  	wooMatrix2D<double> testMatrix(3, 4, simpleData);
  	
  	cout << "Original matrix:" << endl;
  	PrintMatrix(testMatrix);
  	cout << endl;
  	
  	wooMatrix2D<double> testMatrixT = testMatrix.transpose();
  	cout << "The transposed result:" << endl;
  	PrintMatrix(testMatrixT);
  	
  	cout << endl;
  	cout << "Verify the dimensions." << endl;
  	cout << "Original (rows,cols) = (" << testMatrix.getNumRows() << ", " << testMatrix.getNumCols() << ")." << endl;
  	cout << "Transposed (rows,cols) = (" << testMatrixT.getNumRows() << ", " << testMatrixT.getNumCols() << ")." << endl;
  } 
  
  cout << endl;
  cout << "**********************************************************************" << endl;
  cout << endl;
  cout << "Testing linear least squares." << endl;
  
  // Linear least squares - Test 1.
  {
  	// Define our X matrix.
  	std::vector<double> Xdata = {1.0, 1.0, 1.0, 2.0, 1.0, 3.0};
  	wooMatrix2D<double> X(3, 2, Xdata);
  	
  	// Define our y vector.
  	std::vector<double> Ydata = {2.0, 4.0, 4.0};
  	wooVector<double> y(Ydata);
  	
  	cout << "Testing with X = " << endl;
  	PrintMatrix(X);
  	cout << endl;
  	cout << "And y = " << endl;
  	PrintVector(y);
  	cout << endl;
  	
  	// Compute the parameters of best fit.
  	wooVector<double> result;
  	int test = wooLSQ<double>(X, y, result);
  	
  	cout << "Giving beta-hat = " << endl;
  	PrintVector(result);
  	cout << endl;
  }
  
  // Linear least squares - Test 2.
  cout << "***************************************" << endl;
  cout << "Test with a larger number of equations." << endl;
  {
  	// setup a random number generator.
	  std::random_device myRandomDevice;
      std::mt19937 myRandomGenerator(myRandomDevice());
	  std::uniform_real_distribution<double> myDistribution(-1.0, 1.0);  	
	  
	  // setup the test data.
	  int numPoints = 100;
	  double m = 1.5;
	  double c = 0.5;
	  double xMax = 10.0;
	  double xStep = xMax / static_cast<double>(numPoints);
	  
	  wooMatrix2D<double> X(numPoints, 2);
	  wooVector<double> y(numPoints);
      int count = 0;
	  for (double x=0.0; x<xMax; x+=xStep)
	  {
	  	double randomNumber = myDistribution(myRandomGenerator);
        //cout << randomNumber << endl;
	  	X.setElement(count, 0, 1.0);
	  	X.setElement(count, 1, x);
	  	y.setElement(count, ((m*x)+c) + randomNumber);
	  	count++;
	  }
	  
	  // Display these.
	  cout << "X = " << endl;
	  PrintMatrix(X);
	  cout << endl;
	  cout << "y = " << endl;
	  PrintVector(y);
	  cout << endl;
	  
	  // Write out to a CSV file.
	  /*std::ofstream testFile;
	  testFile.open("testFile.csv");
	  for (int i=0; i<numPoints; ++i)
	  {
	  	testFile << X.getElement(i, 1) << "," << y.getElement(i) << std::endl;
	  }
	  testFile.close();*/
	  
	  //  Apply the linear least squares method.
	  wooVector<double> betaHat;
	  int test = wooLSQ<double>(X, y, betaHat);
	  
	  cout << "Giving beta-hat = " << endl;
	  PrintVector(betaHat);
	  cout << endl;
	  
  }
  
	return 0;
}   