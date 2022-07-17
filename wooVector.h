#ifndef WOOVECTOR_H
#define WOOVECTOR_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

template <class T>
class wooVector
{
    public:
    //constructors
    wooVector();
    wooVector(std::vector<T> input);

    //destructor
    ~wooVector();

    //accessor
    int getNumDims();
    T getElement(int index);

    //overloaded operators
    wooVector<T> operator+ (const wooVector<T> &rhs) const;
    wooVector<T> operator- (const wooVector<T> &rhs) const;
    wooVector<T> operator* (const T &rhs) const;

    template <class U> friend wooVector<U> operator* (const U &lhs, const wooVector<U> &rhs);

    //static functions for dot and cross product
    static T dot(const wooVector<T> &a, const wooVector<T> &b);
    static wooVector<T> cross(const wooVector<T> &a, const wooVector<T> &b);

    private:
    std::vector<T> m_vectorData;
    int m_nDims;
};

/************************************************
    CONSTRUCTORS / DESTRUCTOR
*************************************************/

template <class T>
wooVector<T>::wooVector()
{
    m_nDims = 0;
    m_vectorData = std::vector<T>();
}

template <class T>
wooVector<T>::wooVector(std::vector<T> input)
{
    m_nDims = input.size();
    m_vectorData = input;
}

template <class T>
wooVector<T>::~wooVector()
{
    //nothing to destroy
}

/************************************************
    ACCESSORS
*************************************************/
template <class T>
int wooVector<T>::getNumDims()
{
    return m_nDims;
}

template <class T>
T wooVector<T>::getElement(int index)
{
    return m_vectorData[index];
}

/************************************************
    OVERLOADED OPERATORS
*************************************************/
template <class T>
wooVector<T> wooVector<T>::operator+ (const wooVector<T> &rhs) const
{
    //confirm matching dimensions
    if (m_nDims != rhs.m_nDims)
    {
        throw std::invalid_argument("Vector dimensions must match");
    }

    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; i++)
    {
        resultData.push_back(m_vectorData.at(i) + rhs.m_vectorData.at(i));
    }

    wooVector<T> result(resultData);
    return result;
}

template <class T>
wooVector<T> wooVector<T>::operator- (const wooVector<T> &rhs) const
{
    //confirm matching dimensions
    if (m_nDims != rhs.m_nDims)
    {
        throw std::invalid_argument("Vector dimensions must match");
    }

    std::vector<T> resultData;
    for (int i = 0; i < m_nDims; i++)
    {
        resultData.push_back(m_vectorData.at(i) - rhs.m_vectorData.at(i));
    }

    wooVector<T> result(resultData);
    return result;
}

/************************************************
    OVERLOADED OPERATORS (FRIEND)
*************************************************/
template <class T>
wooVector<T> operator* (const T &lhs, const wooVector<T> &rhs)
{
    std::vector<T> resultData;
    for (int i = 0; i < rhs.m_nDims; i++)
    {
        resultData.push_back(lhs * rhs.m_vectorData.at(i));
    }

    wooVector<T> result(resultData);
    return result;
}

/************************************************
    OVERLOADED OPERATORS (FRIEND)
*************************************************/
template <class T>
T wooVector<T>::dot (const wooVector<T> &a, const wooVector<T> &b)
{
    //confirm matching dimensions
    if (a.m_nDims != b.m_nDims)
    {
        throw std::invalid_argument("Vector dimensions must match");
    }

    T cumulativeSum = static_cast<T>(0.0);
    for (int i = 0; i < a.m_nDims; i++)
    {
        cumulativeSum += a.m_vectorData.at(i) * b.m_vectorData.at(i);
    }

    return cumulativeSum;
}

//cross product (only defined here for 3 dimensions)
template <class T>
wooVector<T> wooVector<T>::cross (const wooVector<T> &a, const wooVector<T> &b)
{
    //confirm matching dimensions
    if (a.m_nDims != b.m_nDims)
    {
        throw std::invalid_argument("Vector dimensions must match");
    }

    if (a.m_nDims != 3)
    {
        throw std::invalid_argument("Vectors must be 3 dimensional");
    }

    std::vector<T> resultData;
    resultData.push_back((a.m_vectorData.at(1) * b.m_vectorData.at(2)) - (a.m_vectorData.at(2) * b.m_vectorData.at(1)));
	resultData.push_back(-((a.m_vectorData.at(0) * b.m_vectorData.at(2)) - (a.m_vectorData.at(2) * b.m_vectorData.at(0))));
	resultData.push_back((a.m_vectorData.at(0) * b.m_vectorData.at(1)) - (a.m_vectorData.at(1) * b.m_vectorData.at(0)));

	wooVector<T> result(resultData);
	return result;
}

#endif