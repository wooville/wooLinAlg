#ifndef WOOLINSOLVE_H
#define WOOLINSOLVE_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "wooMatrix2D.h"
#include "wooVector.h"

constexpr int WOOLINSOLVE_NOUNIQUESOLUTION = -1;
constexpr int WOOLINSOLVE_NOSOLUTIONS = -2;

template <typename T>
int wooLinSolve(const wooMatrix2D<T> &aMatrix, const wooVector<T> &bVector, wooVector<T> &resultVector)
{
    wooMatrix2D<T> inputMatrix = aMatrix;

    int origRank = inputMatrix.rank();

    //combine inputs for gaussian elimination
    //extract bvectordata, create matrix with it and join with copy of amatrix
    int numDims = bVector.getNumDims();
    std::vector<T> bVectorData;
    for (int i = 0; i < numDims; i++)
    {
        bVectorData.push_back(bVector.getElement(i));
    }

    wooMatrix2D<T> bMatrix(numDims, 1, &bVectorData);

    inputMatrix.join(bMatrix);

    //gaussian elimination
    wooMatrix2D<T> rowEchMatrix = inputMatrix.rowEchelon();
    
    //avoid redundant gaussian elim by taking rank of row ech'ed matrix
    int augRank = rowEchMatrix.rank();

    //check unique/infinite/no solutions with augRank
    if ((origRank == augRank) && (origRank < inputMatrix.getNumRows()))
    {
        return WOOLINSOLVE_NOUNIQUESOLUTION;
    }
    else if (origRank < augRank)
    {
        return WOOLINSOLVE_NOSOLUTIONS;
    }
    else
    {
        //use back substitution to compute result
        wooVector<T> output(bVectorData);

        int numRows = rowEchMatrix.getNumRows();
        int numCols = rowEchMatrix.getNumCols();

        //loop through rows in reverse
        T currentResult;
        T cumulativeSum;
        T ans;
        int startRow = numRows-1;
        for (int i = startRow; i  >= 0; i--)
        {
            currentResult = rowEchMatrix.getElement(i, numCols-1);

            //cumulativeSum will remain 0 for first loop (to get our first substitution)
            cumulativeSum = static_cast<T>(0.0);
            for (int j = i+1; j < numRows; j++)
            {
                cumulativeSum += (rowEchMatrix.getElement(i,j) * output.getElement(j));
            }

            ans = (currentResult - cumulativeSum) / rowEchMatrix.getElement(i,i);
            output.setElement(i, ans);
        }

        resultVector = output;
    }
    return 1;
}

#endif