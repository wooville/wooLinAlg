/* *************************************************************************************************
	qbMatrixTest

	Code to test the basic functionality of the qbMatrix class contained in qbMatrix.h

	*** INPUTS ***

	None

	*** OUTPUTS ***

	INT				Flag indicating success or failure of the process.
	Created as part of the qbLinAlg linear algebra library, which is intended to be primarily for
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
	cout << "Code to test wooMatrix2D" << endl;

	// *******************************************************************
	// Create an instance of the wooMatrix2D class.
	// This will contain a simple 2D 3x4 matrix
	double simpleData[12] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
	wooMatrix2D<double> testMatrix(3, 4, simpleData);

	// Extract and print the elements of testMatrix.
	cout << endl << "**************************" << endl;
	cout << "3x4 matrix test (testMatrix)." << endl;
	printMatrix(testMatrix);

	// *******************************************************************
	// Test element retrieval.
	cout << endl << "**************************" << endl;
	cout << "Test element retrieval." << endl;
	cout << "Element (0,0) = " << testMatrix.getElement(0, 0) << endl;
	cout << "Element (1,0) = " << testMatrix.getElement(1, 0) << endl;
	cout << "Element (2,0) = " << testMatrix.getElement(2, 0) << endl;
	cout << "Element (0,1) = " << testMatrix.getElement(0, 1) << endl;
	cout << "Element (1,1) = " << testMatrix.getElement(1, 1) << endl;
	cout << "Element (2,1) = " << testMatrix.getElement(2, 1) << endl;
	cout << "Element (5,5) = " << testMatrix.getElement(5, 5) << endl;

	// *******************************************************************
	// Test matrix multiplication.
	cout << endl << "**************************" << endl;
	cout << "Test matrix multiplication." << endl;
	double simpleData2[12] = { 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0 };
	wooMatrix2D<double> testMatrix2(4, 3, simpleData2);
	cout << "4x3 matrix (testMatrix2)" << endl;
	printMatrix(testMatrix2);
	cout << "Multiplication (testMatrix * testMatrix2) result:" << endl;
	wooMatrix2D<double> multTest1 = testMatrix * testMatrix2;
	printMatrix(multTest1);

	cout << endl << "**************************" << endl;
	cout << "testMatrix2 * testMatrix:" << endl;
	wooMatrix2D<double> multTest2 = testMatrix2 * testMatrix;
	printMatrix(multTest2);

	cout << endl << "**************************" << endl;
	cout << "Test multiplication of column vector by matrix." << endl;
	double columnData[3] = { 1.5, 2.5, 3.5 };
	double squareData[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
	wooMatrix2D<double> testColumn(3, 1, columnData);
	wooMatrix2D<double> squareMatrix(3, 3, squareData);
	cout << "Column vector = " << endl;
	printMatrix(testColumn);
	cout << "Square matrix = " << endl;
	printMatrix(squareMatrix);
	cout << "Column vector * Square matrix = " << endl;
	printMatrix(testColumn * squareMatrix);
	cout << "Square matrix * Column vector = " << endl;
	printMatrix(squareMatrix * testColumn);
	cout << "Square matrix + 1.0 = " << endl;
	wooMatrix2D<double> squareMatrix2 = squareMatrix + 1.0;
	printMatrix(squareMatrix2);
	cout << "(Square matrix + 1.0) * Column vector = " << endl;
	printMatrix(squareMatrix2 * testColumn);

	// *******************************************************************
	// Test equality operator
	cout << endl << "**************************" << endl;
	cout << "Test equility operator." << endl;
	cout << "testMatrix == testMatrix2: " << (testMatrix == testMatrix2) << endl;
	cout << "testMatrix2 == testMatrix: " << (testMatrix2 == testMatrix) << endl;
	cout << "(Let testMatrix3 = testMatrix)" << endl;
	wooMatrix2D<double> testMatrix3 = testMatrix;
	cout << "testMatrix == testMatrix3: " << (testMatrix == testMatrix3) << endl;
	cout << "testMatrix3 == testMatrix: " << (testMatrix3 == testMatrix) << endl;

	// *******************************************************************
	// Test matrix addition by scaler.
	cout << endl << "**************************" << endl;
	cout << "Test addition by scaler" << endl;
	cout << "testMatrix + 2.0 = " << endl;
	printMatrix(testMatrix + 2.0);
	cout << endl;
	cout << "2.0 + testMatrix = " << endl;
	printMatrix(2.0 + testMatrix);

	// *******************************************************************
	// Test matrix subtraction by scaler.
	cout << endl << "**************************" << endl;
	cout << "Test addition by scaler" << endl;
	cout << "testMatrix - 2.0 = " << endl;
	printMatrix(testMatrix - 2.0);
	cout << endl;
	cout << "2.0 - testMatrix = " << endl;
	printMatrix(2.0 - testMatrix);

	// *******************************************************************
	// Test matrix multiplication by scaler.
	cout << endl << "**************************" << endl;
	cout << "Test multiplication by scaler" << endl;
	cout << "testMatrix * 2.0 = " << endl;
	printMatrix(testMatrix * 2.0);
	cout << endl;
	cout << "2.0 * testMatrix = " << endl;
	printMatrix(2.0 * testMatrix);

	// *******************************************************************
	// Test formation of identity matrix.
	cout << endl << "**************************" << endl;
	cout << "Test formation of identity matrix." << endl;
	wooMatrix2D<double> identityTest(5, 5);
	identityTest.setToIdentity();
	printMatrix(identityTest);

	// *******************************************************************
	// Test joining of two matrices.
	cout << endl << "**************************" << endl;
	cout << "Test joining of two matrices." << endl;
	wooMatrix2D<double> bigSquare(5, 5);
	bigSquare.join(identityTest);
	printMatrix(bigSquare);

	// *******************************************************************
	// Test matrix inversion
	cout << endl << "**************************" << endl;
	cout << "Test matrix inversion" << endl;
	cout << "Attempt to invert a non-square matrix:" << endl;
	try
	{
		testMatrix.inverse();
	}
	catch (invalid_argument& e)
	{
		cerr << e.what() << endl;
		cout << "Execution would normally stop here with a return code of -1." << endl;
	}

	cout << endl << "**************************" << endl;
	cout << "Test matrix inversion" << endl;
	cout << "Attempt to invert a square matrix:" << endl;
	//double invertTestData[9] = {2.0, 1.0, 1.0, 1.0, 2.0, 3.0, 0.0, 3.0, 1.0};
	double invertTestData[9] = { 2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0 };
	wooMatrix2D<double> invertTest(3, 3, invertTestData);
	wooMatrix2D<double> invertResult = invertTest;
	invertResult.inverse();
	cout << "From:" << endl;
	printMatrix(invertTest);
	cout << "To:" << endl;
	printMatrix(invertResult);

	// *******************************************************************
	// Test inversion of a bigger matrix.
	cout << endl << "**************************" << endl;
	cout << "Test inversion of a bigger matrix." << endl;
	double invertTestData2[25] =
	{ 2.0, 3.0, 4.0, 5.0, 6.0,
	 1.0, 2.0, 3.0, 4.0, 5.0,
	 9.0, 5.0, 3.0, 2.0, 6.0,
	 2.0, 4.0, 6.0, 5.0, 1.0,
	 1.0, 7.0, 5.0, 2.0, 3.0 };
	wooMatrix2D<double> invertTest2(5, 5, invertTestData2);
	wooMatrix2D<double> invertResult2 = invertTest2;
	invertResult2.inverse();
	cout << "From:" << endl;
	printMatrix(invertTest2);
	cout << "To:" << endl;
	printMatrix(invertResult2);

	cout << endl;

	cout << endl << "**************************" << endl;
	cout << "Test multiplication of matrix by it's inverse." << endl;
	cout << "Using invertTest2 * invertResult2:" << endl;
	wooMatrix2D<double> invertAccuracy3 = invertTest2 * invertResult2;
	printMatrix(invertAccuracy3);

	cout << endl;
	cout << "Using invertTest * invertResult:" << endl;
	wooMatrix2D<double> invertAccuracy = invertTest * invertResult;
	printMatrix(invertAccuracy);

	return 0;
}