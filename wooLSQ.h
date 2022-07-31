#ifndef WOOLSQ_H
#define WOOLSQ_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "wooVector.h"
#include "wooMatrix2D.h"

constexpr int WOOLSQ_NOINVERSE = -1;

template <typename T>
int wooLSQ(const wooMatrix2D<T> &Xin, const wooVector<T> &yin, wooVector<T> &result)
{
    //copy X and y
    wooMatrix2D<T> X = Xin;
    wooVector<T> y = yin;

    wooMatrix2D<T> XT = X.transpose();
    wooMatrix2D<T> XTX = XT * X;

    //compute inverse in place or return error
    if (!XTX.inverse())
    {
        return WOOLSQ_NOINVERSE;
    }

    //multiply inverse by XT and y for final result
    wooMatrix2D<T> XTXXT = (XTX * XTz);
    result = (XTXXT * y);

    return 1;
}

#endif