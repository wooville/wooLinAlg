#ifndef MATRIX2D_H
#define MATRIX2D_H

/* *************************************************************************************************
	Matrix2D

	Class to handle two-dimensional matrices and their operations.
	Part of the LinAlg linear algebra library, which was created following the excellent
	and informative series on the YouTube channel QuantitativeBytes:

	www.youtube.com/c/QuantitativeBytes
************************************************************************************************* */

//template to allow matrices of different (presumed numeric) data types
//using Template means that we need to define as well as declare everything in the header file (hence no Matrix2D.cpp)
template <class T>

class Matrix2D
{
public:
	//constructors
	Matrix2D();
	Matrix2D(int nRows, int nCols);
    Matrix2D(int nRows, int nCols, const T *inputData);
	Matrix2D(const Matrix2D<T>& inputMatrix);
    Matrix2D(int nRows, int nCols, const std::vector<T>* inputData);

    //destructor
    ~Matrix2D();

    //config
    bool resize(int numRows, int numCols);
    void setToIdentity();

    //access
    T getElement(int row, int col);
    bool setElement(int row, int col, T elementValue);
    int getNumRows();
    int getNumCols();

    //manipulation
    bool inverse();
    T determinant();

    // Overload == operator.
    bool operator== (const Matrix2D<T>& rhs);
    bool compare (const Matrix2D<T>& matrix1, double tolerance);

    // Overload the assignment operator.
    Matrix2D<T> operator= (const Matrix2D<T>& rhs);

    // Overload +, - and * operators
    // friends to grant access to Matrix2D instance members/functions, since these are not part of the class
    template <class U> friend Matrix2D<U> operator+ (const Matrix2D<U>& lhs, const Matrix2D<U>& rhs);
    template <class U> friend Matrix2D<U> operator+ (const U& lhs, const Matrix2D<U>& rhs);
    template <class U> friend Matrix2D<U> operator+ (const Matrix2D<U>& lhs, const U& rhs);

    template <class U> friend Matrix2D<U> operator- (const Matrix2D<U>& lhs, const Matrix2D<U>& rhs);
    template <class U> friend Matrix2D<U> operator- (const U& lhs, const Matrix2D<U>& rhs);
    template <class U> friend Matrix2D<U> operator- (const Matrix2D<U>& lhs, const U& rhs);

    template <class U> friend Matrix2D<U> operator* (const Matrix2D<U>& lhs, const Matrix2D<U>& rhs);
    template <class U> friend Matrix2D<U> operator* (const U& lhs, const Matrix2D<U>& rhs);
    template <class U> friend Matrix2D<U> operator* (const Matrix2D<U>& lhs, const U& rhs);

    bool separate(Matrix2D<T> *matrix1, Matrix2D<T> *matrix2, int colNum);
    
//private:
public:
    int sub2Ind(int row, int col);
    bool isSquare();
    bool closeEnough(T f1, T f2);
    void swapRow(int i, int j);
    void multAdd(int i, int j, T multFactor);
    void multRow(int i, T multFactor);
    bool join(const Matrix2D<T>& matrix2);
    int findRowWithMaxElement(int colNum, int startingRow);
    Matrix2D<T> findSubMatrix(int colNum, int rowNum);
    void printMatrix();

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
Matrix2D<T>::Matrix2D()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

//constructor given rows/cols (create empty matrix)
template <class T>
Matrix2D<T>::Matrix2D(int nRows, int nCols)
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
Matrix2D<T>::Matrix2D(int nRows, int nCols, const T *inputData)
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
Matrix2D<T>::Matrix2D(int nRows, int nCols, const std::vector<T> *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
        m_matrixData[i] = inputData->at(i);
}

//copy constructor
template <class T>
Matrix2D<T>::Matrix2D(const Matrix2D<T>& inputMatrix)
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
Matrix2D<T>::~Matrix2D()
{
    if (m_matrixData)
        delete[] m_matrixData;

    m_matrixData = nullptr;
}

/************************************************
    CONFIGURATION
*************************************************/
template <class T>
bool Matrix2D<T>::resize(int nCols, int nRows)
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
void Matrix2D<T>::setToIdentity()
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
T Matrix2D<T>::getElement(int row, int col)
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
bool Matrix2D<T>::setElement(int row, int col, T value)
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
int Matrix2D<T>::getNumRows()
{
    return m_nRows;
}

template<class T>
int Matrix2D<T>::getNumCols()
{
    return m_nCols;
}

//inverse using gauss-jordan elimination
template<class T>
bool Matrix2D<T>::inverse()
{
    //can only compute inverse of square matrix
    if (!isSquare())
    {
        throw std::invalid_argument("Matrix must be square to compute inverse");
    }

    //form identity matrix with same dimensions & join it
    Matrix2D<T> identityMatrix(m_nRows, m_nCols);
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
                std::cout << "Swap rows " << cRow << " and " << maxIndex << std::endl;
                swapRow(cRow, maxIndex);
            }

            //set value at (cRow, cCol) = 1 if it isnt already
            if (m_matrixData[sub2Ind(cRow, cCol)] != 1.0)
            {
                T multFactor = 1.0 / m_matrixData[sub2Ind(cRow, cCol)];
                multRow(cRow, multFactor);
                std::cout << "Multiply row " << cRow << " by " << multFactor << std::endl;
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

                    if (!closeEnough(rowOneValue, 0.0))
                    {
                        T correctionFactor = -(curElementValue / rowOneValue);
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                        std::cout << "Multiply row " << rowOneIndex << " by " << correctionFactor <<
                            " and add to row " << rowIndex << std::endl;
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
                        std::cout << "Multiply row " << rowOneIndex << " by " << correctionFactor <<
                            " and add to row " << cRow << std::endl;
                    }
                }
            }
        }

        //separate into left and right halves
        Matrix2D<T> leftHalf;
        Matrix2D<T> rightHalf;
        this->separate(&leftHalf, &rightHalf, originalNumCols);

        //after separating, left half should be identity matrix
        if (leftHalf == identityMatrix)
        {
            completeFlag = true;

            //rebuild matrix with just the right half (inverse)
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

    return completeFlag;
}

//recursive function to calculate matrix determinant
template<class T>
inline T Matrix2D<T>::determinant()
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
            Matrix2D<T> subMatrix = findSubMatrix(0,j);

            cumulativeSum += this->getElement(0, j) * subMatrix.determinant() * sign;
            sign = -sign;
        }
        determinant = cumulativeSum;
    }

    return determinant;
}

/************************************************
    OVERLOADED OPERATORS
*************************************************/

/************************************************
    + OPERATOR
*************************************************/
//matrix + matrix
template <class T>
Matrix2D<T> operator+ (const Matrix2D<T>& lhs, const Matrix2D<T>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows*numCols;
    T* tmpResult = new T[numElements];

    for (int i=0; i<numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//scalar + matrix
template <class T>
Matrix2D<T> operator+ (const T& lhs, const Matrix2D<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs + rhs.m_matrixData[i];
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix + scalar
template <class T>
Matrix2D<T> operator+ (const Matrix2D<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] + rhs;
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

/************************************************
    - OPERATOR
*************************************************/
//matrix - matrix
template <class T>
Matrix2D<T> operator- (const Matrix2D<T>& lhs, const Matrix2D<T>& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//scalar - matrix
template <class T>
Matrix2D<T> operator- (const T& lhs, const Matrix2D<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs - rhs.m_matrixData[i];
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix - scalar
template <class T>
Matrix2D<T> operator- (const Matrix2D<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] - rhs;
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

/************************************************
    * OPERATOR
*************************************************/
//matrix * matrix
template <class T>
Matrix2D<T> operator* (const Matrix2D<T>& lhs, const Matrix2D<T>& rhs)
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
        Matrix2D<T> result(l_numRows, r_numCols, tmpResult);
        delete[] tmpResult;
        return result;
    }
    else
    {
        Matrix2D<T> result(1, 1);
        return result;
    }
}

//scalar * matrix
template <class T>
Matrix2D<T> operator* (const T& lhs, const Matrix2D<T>& rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs * rhs.m_matrixData[i];
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

//matrix * scalar
template <class T>
Matrix2D<T> operator* (const Matrix2D<T>& lhs, const T& rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T* tmpResult = new T[numElements];

    for (int i = 0; i < numElements; i++)
    {
        tmpResult[i] = lhs.m_matrixData[i] * rhs;
    }

    Matrix2D<T> result(numRows, numCols, tmpResult);
    delete[] tmpResult;
    return result;
}

/************************************************
    == OPERATOR
*************************************************/
template <class T>
bool Matrix2D<T>::operator== (const Matrix2D<T>& rhs)
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
bool Matrix2D<T>::compare(const Matrix2D<T>& matrix1, double tolerance)
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
bool Matrix2D<T>::separate(Matrix2D<T>* matrix1, Matrix2D<T>* matrix2, int colNum)
{
    //sizes of new matrices
    int numRows = m_nRows;
    int numCols1 = colNum;
    int numCols2 = m_nCols - colNum;

    //resize will also empty matrix
    matrix1->resize(numRows, numCols1);
    matrix2->resize(numRows, numCols1);

    //loop over original and store 
    for (int row = 0; row < m_nRows; row++)
    {
        for (int col = 0; col < m_nCols; col++)
        {
            if (col < colNum)
            {
                matrix1->setElement(row, col, this->getElement(row, col));
            }
            else
            {
                matrix2->setElement(row, col-colNum, this->getElement(row, col));
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
int Matrix2D<T>::sub2Ind(int row, int col)
{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
        return (row * m_nCols) + col;
    else
        return -1;  //outside of index
}

template<class T>
bool Matrix2D<T>::isSquare()
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

template<class T>
bool Matrix2D<T>::closeEnough(T f1, T f2)
{
    return fabs(f1 - f2) < 1e-9;
}

template<class T>
void Matrix2D<T>::swapRow(int i, int j)
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
void Matrix2D<T>::multAdd(int i, int j, T multFactor)
{
    for (int k = 0; k < m_nCols; k++)
    {
        m_matrixData[sub2Ind(i, k)] += (m_matrixData[sub2Ind(j,k)] * multFactor);
    }
}

//multiply row by scalar
template<class T>
void Matrix2D<T>::multRow(int i, T multFactor)
{
    for (int k = 0; k < m_nCols; k++)
    {
        m_matrixData[sub2Ind(i, k)] *= multFactor;
    }
}

//join another matrix to this one
template<class T>
bool Matrix2D<T>::join(const Matrix2D<T>& matrix2)
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
int Matrix2D<T>::findRowWithMaxElement(int colNum, int startingRow)
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
Matrix2D<T> Matrix2D<T>::findSubMatrix(int rowNum, int colNum)
{
    Matrix2D<T> subMatrix(m_nRows-1, m_nCols-1);
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
void Matrix2D<T>::printMatrix()
{
    int nRows = this->getNumRows;
    int nCols = this->getNumCols;
    for (int row = 0; row<nRows; row++)
    {
        for (int col = 0; col < nCols; col++)
        {
            std::cout << std::fixed << std::setprecision(3) << this->getElement(row, col) << "  ";
        }
        std::cout << std::endl;
    }
}

#endif