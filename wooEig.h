#ifndef WOOEIG_H
#define WOOEIG_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "wooVector.h"
#include "wooMatrix2D.h"
#include "wooQR.h"

constexpr int WOOEIG_MATRIXNOTSQUARE = -1;
constexpr int WOOEIG_MATRIXNOTSYMMETRIC = -2;
constexpr int WOOEIG_MAXITERATIONSEXCEEDED = -3;

//estimate real eigenvalues using QR decomp
//only valid for matrices with all real eigenvalues (aka only symmetric matrices)
template <class T>
int wooEIG_QR(const wooMatrix2D<T> &inputMatrix, std::vector<T> &eigenvalues)
{
    wooMatrix2D<T> A = inputMatrix;

    if (!A.isSquare())
    {
        return WOOEIG_MATRIXNOTSQUARE;
    }

    if (!A.isSymmetric())
    {
        return WOOEIG_MATRIXNOTSYMMETRIC;
    }

    int numRows = A.getNumRows();

    wooMatrix2D<T> identityMatrix(numRows, numRows);
    identityMatrix.setToIdentity();

    wooMatrix2D<T> Q (numRows, numRows);
    wooMatrix2D<T> R (numRows, numRows);

    int iterations = 1000;
    int count = 0;
    bool continueFlag = true;
    while ((count < iterations) && continueFlag)
    {
        //compute QR decomp of A
        int returnCode = wooQR<T>(A, Q, R);

        //find next value of A
        A = R * Q;

        //check if A close enough to upper triangular
        if (A.isRowEchelon())
        {
            continueFlag = false;
        }
        count++;
    }

    for (int i = 0; i < numRows; i++)
    {
        eigenvalues.push_back(A.getElement(i,i));
    }


    if (count == iterations)
    {
        return WOOEIG_MAXITERATIONSEXCEEDED;
    }
    else
    {
        return 0;
    }
}

//inverse power iteration method
template <class T>
int wooEIG_InvPIt(const wooMatrix2D<T> &inputMatrix, const T &eigenvalue, wooVector<T> &eigenVector)
{
    wooMatrix2D<T> A = inputMatrix;

    if (!A.isSquare())
    {
        return WOOEIG_MATRIXNOTSQUARE;
    }

    //random generator
    std::random_device myRandomDevice;
    std::mt19937 myRandomGenerator(myRandomDevice());
    std::uniform_int_distribution<int> myDistribution(1.0, 10.0);

    int numRows = A.getNumRows();

    wooMatrix2D<T> identityMatrix(numRows, numRows);
    identityMatrix.setToIdentity();

    //begin computation
    //generate random vector as starting point for approximation
    wooVector<T> v(numRows);
    for (int i = 0; i < numRows; i++)
    {
        v.setElement(i, static_cast<T>(myDistribution(myRandomGenerator)));
    }

    int iterations = 100;
    int count = 0;
    T deltaThreshold = static_cast<T>(1e-9);
    T delta = static_cast<T>(1e6);
    wooVector<T> prevVector(numRows);
    wooMatrix2D<T> tmpMatrix(numRows, numRows);
    
    while ((count < iterations) && (delta > deltaThreshold))
    {
        //for use computing delta
        prevVector = v;

        //find next value of v
        tmpMatrix = A - (eigenvalue * identityMatrix);
        //tmpMatrix.printMatrix();
        //std::cout << std::endl;
        tmpMatrix.inverseFromAdj();
        //std::cout << tmpMatrix.getNumCols() << " || rows " << v.getNumDims() << std::endl;
        v = tmpMatrix * v;
        v.normalize();
        

        //find new delta
        delta = (v - prevVector).norm();
        
        count++;
    }

    eigenVector = v;

    if (count == iterations)
    {
        return WOOEIG_MAXITERATIONSEXCEEDED;
    }
    else
    {
        return 0;
    }
}

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