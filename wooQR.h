#ifndef WOOQR_H
#define WOOQR_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "wooVector.h"
#include "wooMatrix2D.h"

constexpr int WOOQR_MATRIXNOTSQUARE = -1;

template <class T>
int wooQR(const wooMatrix2D<T> &A, wooMatrix2D<T> &Q, wooMatrix2D<T> &R)
{
    wooMatrix2D<T> inputMatrix = A;

    if (!inputMatrix.isSquare())
    {
        return WOOQR_MATRIXNOTSQUARE;
    }

    int numCols = inputMatrix.getNumCols();

    std::vector<wooMatrix2D<T>> Plist;

    for (int j = 0; j < (numCols-1); j++)
    {
        //a1 is column vector for A
        //b1 is the vector we want to reflect a1 onto
        wooVector<T> a1(numCols-j);
        wooVector<T> b1(numCols-j);
        for (int i = j; i < numCols; i++)
        {
            a1.setElement(i-j, inputMatrix.getElement(i,j));
            b1.setElement(i-j, static_cast<T>(0.0));
        }
        b1.setElement(0, static_cast<T>(1.0));

        T a1norm = a1.norm();


        int sgn = -1;
        if (a1.getElement(0) < static_cast<T>(0.0))
        {
            sgn = 1;
        }

        wooVector<T> u = a1 - (sgn * a1norm * b1);
        wooVector<T> n = u.normalized();

        //convert n to matrix and transpose it
        wooMatrix2D<T> nMat(numCols-j, 1);
        for (int i = 0; i < (numCols-j); i++)
        {
            nMat.setElement(i, 0, n.getElement(i));
        }
        wooMatrix2D<T> nMatT = nMat.transpose();
        wooMatrix2D<T> I (numCols-j, numCols-j);
        I.setToIdentity();

        wooMatrix2D<T> Ptemp = I - static_cast<T>(2.0) * nMat * nMatT;

        //form P matrix with original dimensions using Ptemp
        wooMatrix2D<T> P (numCols, numCols);
        P.setToIdentity();
        for (int row=j; row < numCols; row++)
        {
            for (int col = j; col < numCols; col++)
            {
                P.setElement(row, col, Ptemp.getElement(row-j, col-j));
            }
        }

        //append this P to Plist
        Plist.push_back(P);

        //transform input matrix for next loop
        inputMatrix = P * inputMatrix;
    }

    //compute Q
    wooMatrix2D<T> Qmat = Plist.at(0);
    for (int i = 1; i < (numCols-1); i++)
    {
        Qmat = Qmat * Plist.at(i).transpose();
    }

    Q = Qmat;

    //compute R
    int numElements = Plist.size();
    wooMatrix2D<T> Rmat = Plist.at(numElements-1);
    for (int i = (numElements-2); i >= 0; i--)
    {
        Rmat = Rmat * Plist.at(i);
    }
    Rmat = Rmat * A;

    R = Rmat;

    return 1;
}


#endif