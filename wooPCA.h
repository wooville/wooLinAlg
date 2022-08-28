#ifndef WOOPCA_H
#define WOOPCA_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include "wooVector.h"
#include "wooMatrix2D.h"
#include "wooEig.h"

constexpr int WOOPCA_MATRIXNOTSQUARE = -1;
constexpr int WOOPCA_MATRIXNOTSYMMETRIC = -2;

namespace wooPCA {

template <typename T>
std::vector<T> computeColMeans(const wooMatrix2D<T> &inputData)
{
    int numRows = inputData.getNumRows();
    int numCols = inputData.getNumCols();

    std::vector<T> output;

    for (int j = 0; j < numCols; j++)
    {
        T cumulativeSum = static_cast<T>(0.0);
        for (int i = 0; i < numRows; i++)
        {
            cumulativeSum += inputData.getElement(i,j);
        }

        output.push_back(cumulativeSum / static_cast<T>(numRows));
    }

    return output;
}

template <typename T>
void subtractColumnMeans(wooMatrix2D<T> &inputData, std::vector<T> &columnMeans)
{
    int numRows = inputData.getNumRows();
    int numCols = inputData.getNumCols();

    for (int i = 0; i < numRows; i++)
    {
        for (int j = 0; j < numCols; j++)
        {
            inputData.setElement(i,j, inputData.getElement(i,j) - columnMeans.at(j));
        }
    }
}

//compute covariance matrix
template <typename T>
wooMatrix2D<T> computeCovariance(const wooMatrix2D<T> &X)
{
    //note: X'X rather than traditional XX' because of how the data is required to be arranged
    //one column p for each variable with one row k for each observation
    int numRows = X.getNumRows();
    wooMatrix2D<T> covX = (static_cast<T>(1.0) / static_cast<T>(numRows-1)) * (X.transpose() * X);
    return covX;
}

template <typename T>
int computeEigenvectors(const wooMatrix2D<T> & covarianceMatrix, wooMatrix2D<T> &eigenvectors)
{
    wooMatrix2D<T> X = covarianceMatrix;

    if (!X.isSquare())
    {
        return WOOPCA_MATRIXNOTSQUARE;
    }

    if (!X.isSymmetric())
    {
        return WOOPCA_MATRIXNOTSYMMETRIC;
    }

    std::vector<T> eigenvalues;
    int returnCode = wooEIG_QR(X, eigenvalues);

    std::sort(eigenvalues.begin(), eigenvalues.end());
    std::reverse(eigenvalues.begin(), eigenvalues.end());

    //compute eigenvector for each eigenvalue
    wooVector<T> eV(X.getNumCols());
    wooMatrix2D<T> eVM(X.getNumRows(), X.getNumCols());
    for (int j = 0; j < eigenvalues.size(); j++)
    {
        T eig = eigenvalues.at(j);
        //std::cout << j;
        int returnCode2 = wooEIG_InvPIt<T>(X, eig, eV);
        for (int i = 0; i < eV.getNumDims(); i++)
        {
            eVM.setElement(i, j, eV.getElement(i));
        }
    }

    eigenvectors = eVM;

    return returnCode;
}

//compute principal components
template <typename T>
int wooPCA(const wooMatrix2D<T> &inputData, wooMatrix2D<T> &outputComponents)
{
    wooMatrix2D<T> X = inputData;

    std::vector<T> columnMeans = computeColMeans(inputData);
    subtractColumnMeans<T>(X, columnMeans);
    wooMatrix2D<T> covX = computeCovariance(X);

    wooMatrix2D<T> eigenvectors;
    int returnCode = computeEigenvectors(covX, eigenvectors);

    return returnCode;
}

}

#endif