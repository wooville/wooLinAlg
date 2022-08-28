#ifndef WOOMATRIX2D_H
#define WOOMATRIX2D_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "wooVector.h"

/* *************************************************************************************************
	wooMatrix2D

	Class to handle two-dimensional matrices and their operations.
	Part of the LinAlg linear algebra library, which was created following the excellent
	and informative series on the YouTube channel QuantitativeBytes:

	www.youtube.com/c/QuantitativeBytes
************************************************************************************************* */

//template to allow matrices of different (presumed numeric) data types
//using Template means that we need to define as well as declare everything in the header file (hence no wooMatrix2D.cpp)
template <class T>

class wooMatrix2D
{
public:
	//constructors
	wooMatrix2D();
	wooMatrix2D(int nRows, int nCols);
    wooMatrix2D(int nRows, int nCols, const T *inputData);
	wooMatrix2D(const wooMatrix2D<T> &inputMatrix);
    wooMatrix2D(int nRows, int nCols, const std::vector<T> &inputData);

    //destructor
    ~wooMatrix2D();

    //config
    bool resize(int numRows, int numCols);
    void setToIdentity();

    //access
    T getElement(int row, int col) const;
    bool setElement(int row, int col, T elementValue);
    int getNumRows() const;
    int getNumCols() const; 

    //manipulation
    bool inverse();
    bool inverseFromAdj();
    T determinant();
    wooMatrix2D<T> rowEchelon();
    wooMatrix2D<T> transpose() const;

    wooMatrix2D<T> getCofactor(int p, int q);
    wooMatrix2D<T> adjoint();

    // Overload == operator.
    bool operator== (const wooMatrix2D<T>& rhs);
    bool compare (const wooMatrix2D<T>& matrix1, double tolerance);

    // Overload the assignment operator.
    wooMatrix2D<T> operator= (const wooMatrix2D<T>& rhs);

    // Overload +, - and * operators
    // friends to grant access to wooMatrix2D instance members/functions, since these are not part of the class
    template <class U> friend wooMatrix2D<U> operator+ (const wooMatrix2D<U>& lhs, const wooMatrix2D<U>& rhs);
    template <class U> friend wooMatrix2D<U> operator+ (const U& lhs, const wooMatrix2D<U>& rhs);
    template <class U> friend wooMatrix2D<U> operator+ (const wooMatrix2D<U>& lhs, const U& rhs);

    template <class U> friend wooMatrix2D<U> operator- (const wooMatrix2D<U>& lhs, const wooMatrix2D<U>& rhs);
    template <class U> friend wooMatrix2D<U> operator- (const U& lhs, const wooMatrix2D<U>& rhs);
    template <class U> friend wooMatrix2D<U> operator- (const wooMatrix2D<U>& lhs, const U& rhs);

    template <class U> friend wooMatrix2D<U> operator* (const wooMatrix2D<U>& lhs, const wooMatrix2D<U>& rhs);
    template <class U> friend wooMatrix2D<U> operator* (const U& lhs, const wooMatrix2D<U>& rhs);
    template <class U> friend wooMatrix2D<U> operator* (const wooMatrix2D<U>& lhs, const U& rhs);

    template <class U> friend wooVector<U> operator* (const wooMatrix2D<U>& lhs, const wooVector<U>& rhs);

    bool separate(wooMatrix2D<T> &matrix1, wooMatrix2D<T> &matrix2, int colNum);

    bool isSquare();
    bool isSymmetric();
    bool isRowEchelon();
    bool isNonZero();
    int rank();
    
//private:
public:
    int sub2Ind(int row, int col) const;
    bool closeEnough(T f1, T f2);
    void swapRow(int i, int j);
    void multAdd(int i, int j, T multFactor);
    void multRow(int i, T multFactor);
    bool join(const wooMatrix2D<T>& matrix2);
    int findRowWithMaxElement(int colNum, int startingRow);
    wooMatrix2D<T> findSubMatrix(int colNum, int rowNum);
    void printMatrix();
    void printMatrix(int precision);


private:
    T* m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};

//definitions for above declarations

/************************************************
    CONSTRUCTORS / DESTRUCTOR
*************************************************/

//default constructor
template <class T>
wooMatrix2D<T>::wooMatrix2D()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

//constructor given rows/cols (create empty matrix)
template <class T>
wooMatrix2D<T>::wooMatrix2D(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = 0.0;
}

//constructor given const linear array (create matrix from values)
template <class T>
wooMatrix2D<T>::wooMatrix2D(int nRows, int nCols, const T *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData[i];
}

//constructor given std::vector
template <class T>
wooMatrix2D<T>::wooMatrix2D(int nRows, int nCols, const std::vector<T> &inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData.at(i);
}

//copy constructor
template <class T>
wooMatrix2D<T>::wooMatrix2D(const wooMatrix2D<T>& inputMatrix)
{
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = inputMatrix.m_nElements;

    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputMatrix.m_matrixData[i];
}

//destructor
template <class T>
wooMatrix2D<T>::~wooMatrix2D()
{
    if (m_matrixData)
        delete[] m_matrixData;

    m_matrixData = nullptr;
}

/************************************************
    CONFIGURATION
*************************************************/
template <class T>
bool wooMatrix2D<T>::resize(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = (m_nRows * m_nCols);
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr)
    {
        for (int i = 0; i < m_nElements; i++)
        {
            m_matrixData[i] = 0.0;
        }

        return true;
    }
    else
    {
        return false;
    }
}
//turn this matrix into identity matrix if square
template<class T>
void wooMatrix2D<T>::setToIdentity()
{
    if (!isSquare())
    {
        throw std::invalid_argument("Identity matrix must be square");
    }

    for (int row = 0; row < m_nRows; row++)
    {
        for (int col = 0; col < m_nCols; col++)
        {
            if (col == row)
            {
                m_matrixData[sub2Ind(col, row)] = 1.0;
            }
            else
            {
                m_matrixData[sub2Ind(col, row)] = 0.0;
            }
        }
    }
}

/************************************************
    ACCESSORS
*************************************************/
template <class T>
T wooMatrix2D<T>::getElement(int row, int col) const
{
    int linearIndex = sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        return m_matrixData[linearIndex];
    }
    else
    {
        return 0.0;
    }
}

template <class T>
bool wooMatrix2D<T>::setElement(int row, int col, T value)
{
    int linearIndex = sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        m_matrixData[linearIndex] = value;
        return true;
    }
    else
    {
        return false;
    }
}

template<class T>
int wooMatrix2D<T>::getNumRows() const
{
    return m_nRows;
}

template<class T>
int wooMatrix2D<T>::getNumCols() const
{
    return m_nCols;
}

template<class T>
wooMatrix2D<T> wooMatrix2D<T>::transpose() const
{
    wooMatrix2D<T> resultMatrix(m_nCols, m_nRows);

    for (int i = 0; i < m_nRows; i++)
    {
        for (int j = 0; j < m_nCols; j++)
        {
            resultMatrix.setElement(j, i, this->getElement(i,j));
        }
    }
    
    return resultMatrix;
}

//convert to row echelon using gaussian elimination
template<class T>
wooMatrix2D<T> wooMatrix2D<T>::rowEchelon()
{
    //can only compute inverse of square matrix
    if (m_nCols < m_nRows)
    {
        throw std::invalid_argument("Matrix must have at least as many columns as rows");
    }

    T *tmpMatrixData;
    tmpMatrixData = new T[m_nRows*m_nCols];
    for (int i = 0; i < m_nRows*m_nCols; i++)
    {
        tmpMatrixData[i] = m_matrixData[i];
    }

    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool completeFlag = false;

    while ((!completeFlag) && (count < maxCount))
    {
        for (int diagIndex = 0; diagIndex < m_nRows; diagIndex++)
        {
            //ensure all diag elements = 1
            cRow = diagIndex;
            cCol = diagIndex;

            //set all elements below diagonal to 0
            for (int rowIndex = cRow + 1; rowIndex < m_nRows; rowIndex++)
            {
                //check if 0 already
                if (!closeEnough(m_matrixData[sub2Ind(rowIndex, cCol)], 0.0))
                {
                    int rowOneIndex = cCol;

                    T curElementValue = m_matrixData[sub2Ind(rowIndex, cCol)];
                    T rowOneValue = m_matrixData[sub2Ind(rowOneIndex, cCol)];

                    if (!closeEnough(rowOneValue, 0.0))
                    {
                        T correctionFactor = -(curElementValue / rowOneValue);
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                    }
                }
            }
        }
        completeFlag = true;

        count++;
    }

    wooMatrix2D<T> outputMatrix(m_nRows, m_nCols, m_matrixData);

    for (int i = 0; i < m_nRows*m_nCols; i++)
    {
        m_matrixData[i] = tmpMatrixData[i];
    }

    delete[] tmpMatrixData;

    return outputMatrix;
}

template <class T>
bool wooMatrix2D<T>::isRowEchelon()
{
    T cumulativeSum = 0.0;
    
    for (int i = 1; i < m_nRows; i++)
    {
        for (int j = 0; j < i; j++)
        {
            cumulativeSum += m_matrixData[sub2Ind(i,j)];
        }
    }

    return closeEnough(cumulativeSum, 0.0);
}

//inverse using gauss-jordan elimination
template<class T>
bool wooMatrix2D<T>::inverse()
{
    //can only compute inverse of square matrix
    if (!isSquare())
    {
        throw std::invalid_argument("Matrix must be square to compute inverse");
    }

    //form identity matrix with same dimensions & join it
    wooMatrix2D<T> identityMatrix(m_nRows, m_nCols);
    identityMatrix.setToIdentity();

    int originalNumCols = m_nCols;
    join(identityMatrix);

    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool completeFlag = false;

    while ((!completeFlag) && (count < maxCount))
    {
        for (int diagIndex = 0; diagIndex < m_nRows; diagIndex++)
        {
            //ensure all diag elements = 1
            cRow = diagIndex;
            cCol = diagIndex;

            int maxIndex = findRowWithMaxElement(cCol, cRow);
            
            //if row with max isnt this row, swap into this row
            if (maxIndex != cRow)
            {
                //std::cout << "Swap rows " << cRow << " and " << maxIndex << std::endl;
                swapRow(cRow, maxIndex);
            }

            //set value at (cRow, cCol) = 1 if it isnt already
            //add check != 0 to avoid nan result
            if (m_matrixData[sub2Ind(cRow, cCol)] != 1.0)
            {
                T multFactor = 1.0 / m_matrixData[sub2Ind(cRow, cCol)];
                multRow(cRow, multFactor);
                //std::cout << "Multiply row " << cRow << " by " << multFactor << std::endl;
                //std::cout << "should be 1 " << m_matrixData[sub2Ind(cRow, cCol)] << std::endl;
            }

            //evaluate the column
            for (int rowIndex = cRow + 1; rowIndex < m_nRows; rowIndex++)
            {
                //check if 0 already
                if (!closeEnough(m_matrixData[sub2Ind(rowIndex, cCol)], 0.0))
                {
                    //obtain element of interest from matrix diagonal
                    //should be valid for any invertable matrix
                    int rowOneIndex = cCol;

                    T curElementValue = m_matrixData[sub2Ind(rowIndex, cCol)];
                    T rowOneValue = m_matrixData[sub2Ind(rowOneIndex, cCol)];

                    T correctionFactorCheck = -(curElementValue / rowOneValue);

                    if (!closeEnough(rowOneValue, 0.0))
                    {
                        T correctionFactor = -(curElementValue / rowOneValue);
                        
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                        
                        //std::cout << "111Multiply row " << rowOneIndex << " by " << correctionFactor <<
                        //    " and add to row " << rowIndex << std::endl;
                    }
                }
            }

            //evaluate the row
            for (int colIndex = cCol + 1; colIndex < originalNumCols; colIndex++)
            {
                //check if 0 already
                if (!closeEnough(m_matrixData[sub2Ind(cRow, colIndex)], 0.0))
                {
                    //obtain element of interest from matrix diagonal
                    //should be valid for any invertable matrix
                    int rowOneIndex = colIndex;

                    T curElementValue = m_matrixData[sub2Ind(cRow, colIndex)];
                    T rowOneValue = m_matrixData[sub2Ind(rowOneIndex, colIndex)];

                    if (!closeEnough(rowOneValue, 0.0))
                    {
                        T correctionFactor = -(curElementValue / rowOneValue);
                        multAdd(cRow, rowOneIndex, correctionFactor);
                        //std::cout << "222Multiply row " << rowOneIndex << " by " << correctionFactor <<
                        //    " and add to row " << cRow << std::endl;
                    }
                }
            }
        }

        //separate into left and right halves
        wooMatrix2D<T> leftHalf;
        wooMatrix2D<T> rightHalf;
        
        this->separate(leftHalf, rightHalf, originalNumCols);
        //std::cout << "lh " << leftHalf.getNumCols() << " || rh " << rightHalf.getNumCols() << std::endl;

        //after separating, left half should be identity matrix
        //leftHalf.printMatrix();
        //std::cout << std::endl;
        if (leftHalf == identityMatrix)
        {
            completeFlag = true;

            //rebuild matrix with just the right half (inverse)
            //std::cout << "orig " << originalNumCols << std::endl;
            m_nCols = originalNumCols;
            m_nElements = m_nRows * m_nCols;
            delete[] m_matrixData;
            m_matrixData = new T[m_nElements];
            for (int i = 0; i < m_nElements; i++)
            {
                m_matrixData[i] = rightHalf.m_matrixData[i];
            }
        }

        count++;
    }

    //std::cout << "inside inverse " << m_nCols << std::endl;

    return completeFlag;
}

template<class T>
bool wooMatrix2D<T>::inverseFromAdj()
{
    if (!isSquare())
    {
        throw std::invalid_argument("Can only find inverse of square matrix");
    }

    T det = this->determinant();
    wooMatrix2D<T> adj = this->adjoint();
    wooMatrix2D<T> inverse = (1/det)*adj;

    //update stored data
    for (int i = 0; i < m_nElements; i++)
    {
        m_matrixData[i] = inverse.m_matrixData[i];
    }

    return 1;
}

//recursive function to calculate matrix determinant
template<class T>
T wooMatrix2D<T>::determinant()
{
    if (!isSquare())
    {
        throw std::invalid_argument("Can only find determinant of square matrix");
    }
    //base case 2x2 matrix
    T determinant;
    if (m_nRows == 2)
    {
        determinant = m_matrixData[0] * m_matrixData[3] - m_matrixData[1] * m_matrixData[2];
    }
    else
    {
        T sign = 1.0;
        T cumulativeSum = 0.0;

        for (int j = 0; j < m_nCols; j++)
        {
            wooMatrix2D<T> subMatrix = findSubMatrix(0,j);

            cumulativeSum += this->getElement(0, j) * subMatrix.determinant() * sign;
            sign = -sign;
        }
        determinant = cumulativeSum;
    }

    return determinant;
}

// Function to get cofactor of A[p][q]
template<class T>
wooMatrix2D<T> wooMatrix2D<T>::getCofactor(int p, int q)
{
    int i = 0, j = 0;
    wooMatrix2D<T> tmp(m_nRows-1, m_nCols-1);
 
    // Looping for each element of the matrix
    for (int row = 0; row < m_nRows; row++) {
        for (int col = 0; col < m_nCols; col++) {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q) {
                tmp.setElement(i, j++, m_matrixData[sub2Ind(row, col)]);
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == m_nCols - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }

    return tmp;
}

// Function to get adjoint of this matrix.
template<class T>
wooMatrix2D<T> wooMatrix2D<T>::adjoint()
{
    if (!isSquare())
    {
        throw std::invalid_argument("Can only find determinant of square matrix");
    }

    wooMatrix2D<T> adj(m_nRows, m_nCols);

    if (m_nRows == 1) {
        adj.setElement(0, 0, 1);
        return adj;
    }
 
    // tmp is used to store cofactors
    int sign = 1;
    wooMatrix2D<T> tmp(m_nRows-1, m_nCols-1);
 
    for (int i = 0; i < m_nRows; i++) {
        for (int j = 0; j < m_nCols; j++) {
            // get cofactor of A[i][j]
            tmp = getCofactor(i, j);
 
            // sign of adj[j][i] positive if sum of row and column indexes is even
            sign = ((i + j) % 2 == 0) ? 1 : -1;
 
            //get transpose of cofactor matrix
            T val = (sign) * (tmp.determinant());
            adj.setElement(j, i, val);
        }
    }

    return adj;
}

template <class T>
int wooMatrix2D<T>::rank()
{
    wooMatrix2D<T> matrixCopy = this->rowEchelon();

    int numNonZeroRows = 0;
    if (!matrixCopy.isRowEchelon())
    {
        std::vector<wooMatrix2D<T>> subMatrixVector;
        subMatrixVector.push_back(*this);

        bool completeFlag = false;
        int subMatrixCount = 0;
        while((subMatrixCount < subMatrixVector.size()) && (!completeFlag))
        {
            wooMatrix2D<T> curMatrix = subMatrixVector[subMatrixCount];
            subMatrixCount++;

            if (curMatrix.isNonZero())
            {
                T curMatrixDet = curMatrix.determinant();

                if (!closeEnough(curMatrixDet, 0.0))
                {
                    completeFlag = true;
                    numNonZeroRows = curMatrix.getNumRows();
                }
                else
                {
                    if ((curMatrix.getNumRows() > 2) && (curMatrix.getNumCols() > 2))
                    {
                        for (int i = 0; i < getNumRows(); i++)
                        {
                            for (int j = 0; j < getNumCols(); j++)
                            {
                                subMatrixVector.push_back(curMatrix.findSubMatrix(i,j));
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        int nRows = matrixCopy.getNumRows();
        int nCols = matrixCopy.getNumCols();

        for (int i=0; i<nRows; i++)
        {
            int colSum = 0;
            for (int j=0; j<nCols; j++)
            {
                if (!closeEnough(matrixCopy.getElement(i,j), 0.0))
                {
                    colSum++;
                }
            }

            if (colSum > 0)
            {
                numNonZeroRows++;
            }

        }
    }

    return numNonZeroRows;
}

/************************************************
    OVERLOADED OPERATORS
*************************************************/
/************************************************
    = OPERATOR
*************************************************/
template <class T>
wooMatrix2D<T> wooMatrix2D<T>::operator= (const wooMatrix2D<T>& rhs)
{
	// Make sure we're not assigning to ourself.
	if (this != &rhs)
	{
		m_nRows = rhs.m_nRows;
		m_nCols = rhs.m_nCols;
		m_nElements = rhs.m_nElements;
		
		if (m_matrixData)
			delete[] m_matrixData;
		
		m_matrixData = new T[m_nElements];
		for (int i=0; i<m_nElements; i++)
			m_matrixData[i] = rhs.m_matrixData[i];	
			
		return *this;
	}
}

/************************************************
    + OPERATOR
*************************************************/
//matrix + matrix
template <class T>
wooMatrix2D<T> operator+ (const wooMatrix2D<T>& lhs, const wooMatrix2D<T>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows*numCols;
    T* tmpResult = new T[numElements];

    for (int i=0; i<numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//scalar + matrix
template <class T>
wooMatrix2D<T> operator+ (const T& lhs, const wooMatrix2D<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs + rhs.m_matrixData[i];
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix + scalar
template <class T>
wooMatrix2D<T> operator+ (const wooMatrix2D<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] + rhs;
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

/************************************************
    - OPERATOR
*************************************************/
//matrix - matrix
template <class T>
wooMatrix2D<T> operator- (const wooMatrix2D<T>& lhs, const wooMatrix2D<T>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//scalar - matrix
template <class T>
wooMatrix2D<T> operator- (const T& lhs, const wooMatrix2D<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs - rhs.m_matrixData[i];
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix - scalar
template <class T>
wooMatrix2D<T> operator- (const wooMatrix2D<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] - rhs;
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

/************************************************
    * OPERATOR
*************************************************/
//matrix * matrix
template <class T>
wooMatrix2D<T> operator* (const wooMatrix2D<T>& lhs, const wooMatrix2D<T>& rhs)
{
    int l_numRows = lhs.m_nRows;
    int l_numCols = lhs.m_nCols;
    int r_numRows = rhs.m_nRows;
    int r_numCols = rhs.m_nCols;
    int resultIndex;

    if (l_numCols == r_numRows) //req condition for matrix multiplication
    {
        T* tmpResult = new T[lhs.m_nRows * rhs.m_nCols];

        for (int lhsRow = 0; lhsRow < l_numRows; lhsRow++) //loop through LHS rows
        {
            for (int rhsCol = 0; rhsCol < r_numCols; rhsCol++) //loop through RHS cols
            {
                T elementResult = 0.0;
                for (int lhsCol = 0; lhsCol < l_numCols; lhsCol++) //loop through elements of LHS row
                {
                    int lhsIndex = (lhsRow * l_numCols) + lhsCol;
                    int rhsIndex = (lhsCol * r_numCols) + rhsCol;   //RHS row = LHS col

                    elementResult += (lhs.m_matrixData[lhsIndex] * rhs.m_matrixData[rhsIndex]);
                }

                resultIndex = (lhsRow * r_numCols) + rhsCol;
                tmpResult[resultIndex] = elementResult;
            }
        }
        wooMatrix2D<T> result(l_numRows, r_numCols, tmpResult);
        delete[] tmpResult;
        return result;
    }
    else
    {
        wooMatrix2D<T> result(1, 1);
        return result;
    }
}

//scalar * matrix
template <class T>
wooMatrix2D<T> operator* (const T& lhs, const wooMatrix2D<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs * rhs.m_matrixData[i];
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix * scalar
template <class T>
wooMatrix2D<T> operator* (const wooMatrix2D<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] * rhs;
    }

    wooMatrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix * vector
template <class T>
wooVector<T> operator* (const wooMatrix2D<T>& lhs, const wooVector<T>& rhs)
{
    if (lhs.m_nCols != rhs.getNumDims())
    {
        throw std::invalid_argument("Number of cols in matrix must match number of rows in vector");
    }

    wooVector<T> result = rhs;

    T cumulativeSum;
    for (int row = 0; row < lhs.m_nRows; row++)
    {
        cumulativeSum = static_cast<T>(0.0);
        for (int col = 0; col < lhs.m_nCols; col++)
        {
            cumulativeSum += (lhs.getElement(row,col) * rhs.getElement(col));
        }
        result.setElement(row, cumulativeSum);
    }
    
    return result;
}

/************************************************
    == OPERATOR
*************************************************/
template <class T>
bool wooMatrix2D<T>::operator== (const wooMatrix2D<T>& rhs)
{
    //check size before looping through
    if ((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
    {
        return false;
    }

    bool flag = true;
    for (int i = 0; i < this->m_nElements; i++)
    {
        if (!closeEnough(this->m_matrixData[i], rhs.m_matrixData[i]))
        {
            flag = false;
            break;
        }
    }
    return flag;
}

template<class T>
bool wooMatrix2D<T>::compare(const wooMatrix2D<T>& matrix1, double tolerance)
{
    int numRows1 = matrix1.m_nRows;
    int numCols1 = matrix1.m_nCols;

    if ((numRows1 != m_nRows) || (numCols1 != m_nCols))
    {
        return false;
    }

    //compute sum of squares difference
    double cumulativeSum = 0.0;
    for (int i = 0; i < m_nElements; i++)
    {
        T element1 = matrix1.m_matrixData[i];
        T element2 = m_matrixData[i];
        cumulativeSum += ((element1 - element2) * (element1 - element2));
    }

    double ssd = sqrt(cumulativeSum / ((numRows1 * numCols1) - 1));

    if (ssd < tolerance)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//separate matrix into 2 matrices around colNum
template<class T>
bool wooMatrix2D<T>::separate(wooMatrix2D<T> &matrix1, wooMatrix2D<T> &matrix2, int colNum)
{
    //sizes of new matrices
    int numRows = m_nRows;
    int numCols1 = colNum;
    int numCols2 = m_nCols - colNum;

    //resize will also empty matrix
    matrix1.resize(numRows, numCols1);
    matrix2.resize(numRows, numCols1);

    //loop over original and store 
    for (int row = 0; row < m_nRows; row++)
    {
        for (int col = 0; col < m_nCols; col++)
        {
            if (col < colNum)
            {
                matrix1.setElement(row, col, this->getElement(row, col));
            }
            else
            {
                matrix2.setElement(row, col-colNum, this->getElement(row, col));
            }
        }
    }

    
    return true;
}

/************************************************
    PRIVATE
*************************************************/
//return linear index of element given row/col subscript
template <class T>
int wooMatrix2D<T>::sub2Ind(int row, int col) const
{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
        return (row * m_nCols) + col;
    else
        return -1;  //outside of index
}

template<class T>
bool wooMatrix2D<T>::isSquare()
{
    if (m_nCols == m_nRows)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//test symmetry
template<class T>
bool wooMatrix2D<T>::isSymmetric()
{
    if (!this->isSquare())
    {
        return false;
    }
    

    T curRowElement = static_cast<T>(0.0);
    T curColElement = static_cast<T>(0.0);
    bool returnFlag = true;
    int diagIndex = 0;
    while ((diagIndex < m_nCols) && returnFlag)
    {
        int rowIndex = diagIndex + 1;
        while ((rowIndex < m_nCols) && returnFlag)
        {
            curRowElement = this->getElement(rowIndex, diagIndex);
            curColElement = this->getElement(diagIndex, rowIndex);

            if (!closeEnough(curRowElement, curColElement))
            {
                returnFlag = false;
            }

            rowIndex++;
        }

        diagIndex++;
    }

    return returnFlag;
}

template<class T>
bool wooMatrix2D<T>::closeEnough(T f1, T f2)
{
    return fabs(f1 - f2) < 1e-9;
}

//test by taking sum across lower triangular and checking == 0
template<class T>
bool wooMatrix2D<T>::isNonZero()
{
    int numNonZero = 0;
    for (int i = 0; i < m_nElements; i++)
    {
        if (!closeEnough(m_matrixData[i], 0.0))
        {
            numNonZero++;
        }
    }

    return (numNonZero != 0);
}

template<class T>
void wooMatrix2D<T>::swapRow(int i, int j)
{
    T* tmpRow = new T[m_nCols];
    for (int k = 0; k < m_nCols; k++)
    {
        tmpRow[k] = m_matrixData[sub2Ind(i,k)];
    }

    for (int k = 0; k < m_nCols; k++)
    {
        m_matrixData[sub2Ind(i, k)] = m_matrixData[sub2Ind(j,k)];
    }

    for (int k = 0; k < m_nCols; k++)
    {
        m_matrixData[sub2Ind(j, k)] = tmpRow[k];
    }

    delete[] tmpRow;
}

//add multiple of row j to row i
template<class T>
void wooMatrix2D<T>::multAdd(int i, int j, T multFactor)
{
    for (int k = 0; k < m_nCols; k++)
    {
        m_matrixData[sub2Ind(i, k)] += (m_matrixData[sub2Ind(j,k)] * multFactor);
    }
}

//multiply row by scalar
template<class T>
void wooMatrix2D<T>::multRow(int i, T multFactor)
{
    for (int k = 0; k < m_nCols; k++)
    {
        m_matrixData[sub2Ind(i, k)] *= multFactor;
    }
}

//join another matrix to this one
template<class T>
bool wooMatrix2D<T>::join(const wooMatrix2D<T>& matrix2)
{
    int numRows1 = m_nRows;
    int numRows2 = matrix2.m_nRows;
    int numCols1 = m_nCols;
    int numCols2 = matrix2.m_nCols;

    if (numRows1 != numRows2)
    {
        throw std::invalid_argument("Matrices must have same number of rows");
    }

    //only col size changes
    T* newMatrixData = new T[numRows1*(numCols1+numCols2)];

    int index, resultIndex;
    for (int i = 0; i < numRows1; ++i)
    {
        for (int j = 0; j < (numCols1+numCols2); j++)
        {
            resultIndex = (i * (numCols1+numCols2)) + j;

            //if j is in lhs matrix
            if (j < numCols1)
            {
                index = (i * numCols1) + j;
                newMatrixData[resultIndex] = m_matrixData[index];
            }
            else //rhs matrix
            {
                index = (i * numCols2) + (j - numCols1);
                newMatrixData[resultIndex] = matrix2.m_matrixData[index];
            }
        }
    }

    //update stored data
    m_nCols = numCols1 + numCols2;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
    {
        m_matrixData[i] = newMatrixData[i];
    }
    delete[] newMatrixData;

    return true;
}

template<class T>
int wooMatrix2D<T>::findRowWithMaxElement(int colNum, int startingRow)
{
    T tmp = m_matrixData[sub2Ind(startingRow, colNum)];
    int rowIndex = startingRow;
    for (int k = startingRow+1; k < m_nRows; k++)
    {
        if (fabs(m_matrixData[sub2Ind(k, colNum)]) > fabs(tmp))
        {
            rowIndex = k;
            tmp = m_matrixData[sub2Ind(k, colNum)];
        }
    }

    return rowIndex;
}

template<class T>
wooMatrix2D<T> wooMatrix2D<T>::findSubMatrix(int rowNum, int colNum)
{
    wooMatrix2D<T> subMatrix(m_nRows-1, m_nCols-1);
    int count = 0;
    for (int i = 0; i < m_nRows; i++)
    {
        for (int j = 0; j < m_nCols; j++)
        {
            if (i != rowNum && j != colNum)
            {
                subMatrix.m_matrixData[count] = this->getElement(i, j);
                count++;
            }
        }
    }
    return subMatrix;
}

template<class T>
void wooMatrix2D<T>::printMatrix()
{
    int nRows = this->getNumRows();
    int nCols = this->getNumCols();
    for (int row = 0; row<nRows; row++)
    {
        for (int col = 0; col < nCols; col++)
        {
            std::cout << std::fixed << std::setprecision(3) << this->getElement(row, col) << "  ";
        }
        std::cout << std::endl;
    }
}

template<class T>
void wooMatrix2D<T>::printMatrix(int precision)
{
    int nRows = this->getNumRows();
    int nCols = this->getNumCols();
    for (int row = 0; row<nRows; row++)
    {
        for (int col = 0; col < nCols; col++)
        {
            std::cout << std::fixed << std::setprecision(precision) << this->getElement(row, col) << "  ";
        }
        std::cout << std::endl;
    }
}

#endif