#ifndef WOOEIG_H
#define WOOEIG_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "wooVector.h"
#include "wooMatrix2D.h"

constexpr int WOOEIG_MATRIXNOTSQUARE = -1;

//power iteration method of finding eigenvector/value
template <class T>
int wooEIG_PIt(const wooMatrix2D<T> X, T &eigenvalue, wooVector<T> &eigenVector)
{
    wooMatrix2D<T> inputMatrix = X;

    if (!inputMatrix.isSquare())
    {
        return WOOEIG_MATRIXNOTSQUARE;
    }

    //random generator
    std::random_device myRandomDevice;
    std::mt19937 myRandomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    int numRows = inputMatrix.getNumRows();

    wooMatrix2D<T> identityMatrix(numRows, numRows);
    identityMatrix.setToIdentity();

    //begin computation
    //generate random vector as starting point for approximation
    wooVector<T> v(numRows);
    for (int i = 0; i < numRows; i++)
    {
        v.setElement(i, static_cast<T>(myDistribution(myRandomGenerator)));
    }

    //approximate eigenvector with PI
    wooVector<T> v1(numRows);
    int iterations = 1000;
    for (int i = 0; i < iterations; i++)
    {
        v1 = inputMatrix * v;
        v1.normalize();
        v = v1;
    }

    eigenVector = v1;

    //compute corresponding eigenvalue
    T cumulativeSum = static_cast<T>(0.0);
    for (int i = 1; i < numRows; i++)
    {
        cumulativeSum += inputMatrix.getElement(0,i) * v1.getElement(i);
    }

    eigenvalue = (cumulativeSum / v1.getElement(0)) + inputMatrix.getElement(0,0);

    return 0;
}

#endif