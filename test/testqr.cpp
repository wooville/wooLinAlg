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
#include "../wooQR.h"

using namespace std;

int main()
{
	cout << "**********************************************" << endl;
	cout << "Testing QR decomposition code." << endl;
	cout << "**********************************************" << endl;
	cout << endl;

	{	
		cout << "Testing with simple 3x3 matrix:" << endl;
	
		std::vector<double> simpleData = {0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
		wooMatrix2D<double> testMatrix(3, 3, simpleData);
		
		testMatrix.printMatrix();
		
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		wooMatrix2D<double> Q (3,3);
		wooMatrix2D<double> R (3,3);
		int status = wooQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;
		cout << "QR = " << endl;
		wooMatrix2D<double> QR = Q*R;
		QR.printMatrix();
		cout << endl;		
	}
	
	{
		cout << "Testing with simple 4x4 matrix:" << endl;
		std::vector<double> simpleData = {1.0, 5.0, 3.0, 4.0, 7.0, 8.0, 2.0, 9.0, 7.0, 3.0, 2.0, 1.0, 9.0, 3.0, 5.0, 7.0};
		wooMatrix2D<double> testMatrix(4, 4, simpleData);
		testMatrix.printMatrix();
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		wooMatrix2D<double> Q (4,4);
		wooMatrix2D<double> R (4,4);
		int status = wooQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;		
		cout << "QR = " << endl;
		wooMatrix2D<double> QR = Q*R;
		QR.printMatrix();
		cout << endl;		
	}
	
	{
		cout << "Testing with simple 5x5 matrix:" << endl;
		std::vector<double> simpleData = {2, 6, 4, 6, 8, 6, 7, 9, 7, 9, 2, 3, 6, 3, 5, 6, 1, 1, 5, 5, 3, 5, 6, 5, 6};
		wooMatrix2D<double> testMatrix(5, 5, simpleData);
		testMatrix.printMatrix();
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		wooMatrix2D<double> Q (5,5);
		wooMatrix2D<double> R (5,5);
		int status = wooQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;
		cout << "QR = " << endl;
		wooMatrix2D<double> QR = Q*R;
		QR.printMatrix();
		cout << endl;
	}	
	
	{
		cout << "Testing with simple 5x5 float matrix:" << endl;
		std::vector<double> simpleData = {8.662634278267483, 2.3440981169711796, 3.414158790068152, 9.819959485632891, 9.812414578216162, 4.8096369839436495, 7.743133259609277, 9.871217856632036, 7.100783013043249, 8.127838524397976, 1.3468248609110365, 1.3120774834063536, 9.607366488550678, 2.852679282078192, 8.087038227451359, 7.556075051454403, 5.80117852857823, 3.550189544341768, 3.7807047754393994, 7.934423413357392, 2.866445996919499, 7.125441061546031, 4.53141730712106, 4.297092147605687, 2.5126585000174146};
		wooMatrix2D<double> testMatrix(5, 5, simpleData);
		testMatrix.printMatrix();
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		wooMatrix2D<double> Q (5,5);
		wooMatrix2D<double> R (5,5);
		int status = wooQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;		
		cout << "QR = " << endl;
		wooMatrix2D<double> QR = Q*R;
		QR.printMatrix();
		cout << endl;		
	}		
	
	{
		cout << "Testing with simple 10x10 matrix:" << endl;
		std::vector<double> simpleData = {8, 2, 1, 3, 9, 2, 1, 9, 8, 9, 9, 1, 4, 3, 8, 8, 7, 9, 4, 2, 2, 4, 2, 8, 7, 2, 8, 4, 5, 6, 9, 3, 7, 1, 8, 6, 7, 5, 8, 8, 9, 7, 9, 3, 9, 1, 7, 1, 9, 5, 6, 4, 3, 6, 6, 1, 4, 5, 5, 7, 7, 6, 9, 9, 5, 8, 1, 1, 7, 9, 6, 2, 1, 8, 2, 3, 8, 7, 7, 6, 2, 8, 4, 8, 1, 7, 8, 4, 5, 5, 3, 2, 6, 5, 2, 6, 9, 7, 3, 7};
		wooMatrix2D<double> testMatrix(10, 10, simpleData);
		testMatrix.printMatrix();
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		wooMatrix2D<double> Q (10,10);
		wooMatrix2D<double> R (10,10);
		int status = wooQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;		
		cout << "QR = " << endl;
		wooMatrix2D<double> QR = Q*R;
		QR.printMatrix();
		cout << endl;		
	}		
	
	{
		cout << "Testing with simple 10x10 <float> matrix:" << endl;
		std::vector<float> simpleData = {8, 2, 1, 3, 9, 2, 1, 9, 8, 9, 9, 1, 4, 3, 8, 8, 7, 9, 4, 2, 2, 4, 2, 8, 7, 2, 8, 4, 5, 6, 9, 3, 7, 1, 8, 6, 7, 5, 8, 8, 9, 7, 9, 3, 9, 1, 7, 1, 9, 5, 6, 4, 3, 6, 6, 1, 4, 5, 5, 7, 7, 6, 9, 9, 5, 8, 1, 1, 7, 9, 6, 2, 1, 8, 2, 3, 8, 7, 7, 6, 2, 8, 4, 8, 1, 7, 8, 4, 5, 5, 3, 2, 6, 5, 2, 6, 9, 7, 3, 7};
		wooMatrix2D<float> testMatrix(10, 10, simpleData);
		testMatrix.printMatrix();
		cout << endl;
		cout << "Computing QR decomposition..." << endl;
		wooMatrix2D<float> Q (10,10);
		wooMatrix2D<float> R (10,10);
		int status = wooQR(testMatrix, Q, R);
		cout << "R = " << endl;
		R.printMatrix();
		cout << endl;
		cout << "Q = " << endl;
		Q.printMatrix();
		cout << endl;		
		cout << "QR = " << endl;
		wooMatrix2D<float> QR = Q*R;
		QR.printMatrix();
		cout << endl;		
	}	
	
	return 0;
}