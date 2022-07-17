#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "../wooMatrix2D.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void printMatrix(wooMatrix2D<T> matrix)
{
	int nRows = matrix.getNumRows();
	int nCols = matrix.getNumCols();
	for (int row = 0; row < nRows; ++row)
	{
		for (int col = 0; col < nCols; ++col)
		{
			cout << std::fixed << std::setprecision(3) << matrix.getElement(row, col) << "  ";
		}
		cout << endl;
	}
}

int main()
{

	cout << "Testing implementation of determinant calculation." << endl;
	cout << endl;

	cout << "Generate a test matrix." << endl;
	double testData[9] = { 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0 };
	wooMatrix2D<double> testMatrix(3, 3, testData);
	printMatrix(testMatrix);
	cout << endl;

	cout << "Extract sub-matrix for element (0,0)" << endl;
	wooMatrix2D<double> minor1 = testMatrix.findSubMatrix(0, 0);
	printMatrix(minor1);
	cout << endl;

	cout << "Extract sub-matrix for element (0,1)" << endl;
	wooMatrix2D<double> minor2 = testMatrix.findSubMatrix(0, 1);
	printMatrix(minor2);
	cout << endl;

	cout << "Extract sub-matrix for element (0,2)" << endl;
	wooMatrix2D<double> minor3 = testMatrix.findSubMatrix(0, 2);
	printMatrix(minor3);
	cout << endl;

	cout << "Extract sub-matrix for element (1,1)" << endl;
	wooMatrix2D<double> minor4 = testMatrix.findSubMatrix(1, 1);
	printMatrix(minor4);
	cout << endl;

	cout << "Test with a larger matrix." << endl;
	double testData2[25] =
	{ 2.0, 3.0, 4.0, 5.0, 6.0,
	 1.0, 2.0, 3.0, 4.0, 5.0,
	 9.0, 5.0, 3.0, 2.0, 6.0,
	 2.0, 4.0, 6.0, 5.0, 1.0,
	 1.0, 7.0, 5.0, 2.0, 3.0 };
	wooMatrix2D<double> testMatrix2(5, 5, testData2);
	printMatrix(testMatrix2);
	cout << endl;

	cout << "Extract sub-matrix for element (0,0)" << endl;
	wooMatrix2D<double> minor5 = testMatrix2.findSubMatrix(0, 0);
	printMatrix(minor5);
	cout << endl;

	cout << "Extract sub-matrix for element (0,1)" << endl;
	wooMatrix2D<double> minor6 = testMatrix2.findSubMatrix(0, 1);
	printMatrix(minor6);
	cout << endl;

	cout << "Extract sub-matrix for element (0,2)" << endl;
	wooMatrix2D<double> minor7 = testMatrix2.findSubMatrix(0, 2);
	printMatrix(minor7);
	cout << endl;

	cout << "Extract sub-matrix for element (1,1)" << endl;
	wooMatrix2D<double> minor8 = testMatrix2.findSubMatrix(1, 1);
	printMatrix(minor8);
	cout << endl;

	cout << "Test determinant of 3x3 matrix:" << endl;
	cout << testMatrix.determinant() << endl;
	cout << endl;

	cout << "Test determinant of 5x5 matrix:" << endl;
	cout << testMatrix2.determinant() << endl;
	cout << endl;

	cout << "Test determinant of a singular matrix:" << endl;
	double testData3[9] =
	{ 1.0, 1.0, 1.0,
	 0.0, 1.0, 0.0,
	 1.0, 0.0, 1.0 };
	wooMatrix2D<double> testMatrix3(3, 3, testData3);
	printMatrix(testMatrix3);
	cout << endl;
	cout << "determinant = " << testMatrix3.determinant() << endl;
	cout << endl;

	return 0;
}