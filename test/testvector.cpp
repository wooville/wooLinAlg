#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>

#include "../wooVector.h"

using namespace std;

// A simple function to print a matrix to stdout.
template <class T>
void printVector(wooVector<T> vector)
{
	int nRows = vector.getNumDims();
	for (int row = 0; row < nRows; ++row)
	{
		cout << std::fixed << std::setprecision(3) << vector.getElement(row) << endl;
	}
}

int main()
{

	cout << "Testing implementation of vector class" << endl;
	cout << endl;

	cout << "Generate a test vector." << endl;
	std::vector<double> testData1 = {1.0, 2.0, 3.0};
	wooVector<double> testVector(testData1);
	printVector(testVector);
	cout << endl;

	cout << "test vector arithmetic" << endl;
    cout << "add" << endl;
	std::vector<double> testData2 = {2.0, 4.0, 6.0};
	wooVector<double> testVector2(testData2);
    wooVector<double> testVector3 = testVector + testVector2;
    printVector(testVector3);

    cout << "subtract" << endl;
	testVector3 = testVector - testVector2;
    printVector(testVector3);

    cout << "scalar multiply by 2.0" << endl;
    testVector3 = 2.0*testVector2;
    printVector(testVector3);
    cout << endl;

    cout << "dot product" << endl;
    double dot = wooVector<double>::dot(testVector, testVector2);
    cout << dot << endl;

    cout << "cross product parallel" << endl;
	testVector3 = wooVector<double>::cross(testVector, testVector2);
    printVector(testVector3);
    std::vector<double> testData4 = {-1.0, 5.0, -3.2};
	wooVector<double> testVector4(testData4);
    cout << "cross product non parallel" << endl;
	testVector3 = wooVector<double>::cross(testVector, testVector4);
    printVector(testVector3);
	cout << endl; 

	cout << "test norm" << endl;
	cout << std::setprecision(3) << testVector.norm() << endl;
	cout << endl; 
	cout << "test normalized" << endl;
	wooVector<double> testVectorNorm = testVector.normalized();
	printVector(testVectorNorm);
	cout << "test normalize" << endl;
	testVector.normalize();
	printVector(testVector);

	return 0;
}