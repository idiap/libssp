/*
 * Copyright 2016 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, August 2016
 */

#include "ssp/window.h"

using namespace std;
using namespace ssp;

/**
 * General raised cosine window
 */
static var raisedCosine(
    int iSize, bool iPeriodic, const initializer_list<float> iCoeff
)
{
    static const float pi = atan(1.0) * 4;
    int d = iPeriodic ? iSize : iSize-1;

    var w(iSize, 0.0f);
    for (int i=0; i<iSize; i++)
    {
        int m = -1;
        int j = 0;
        for (const float* c=begin(iCoeff); c!=end(iCoeff); ++c)
        {
            m *= -1;
            w[i] += *c * cos(2.0 * pi * i * j++ / d) * m;
        }
    }
    return w;
}

void Hann::set(int iSize, bool iPeriodic, float iParam)
{
    mWindow = raisedCosine(iSize, iPeriodic, {0.5, 0.5});
}

void Hamming::set(int iSize, bool iPeriodic, float iParam)
{
    mWindow = raisedCosine(iSize, iPeriodic, {0.54, 0.46});
}

void Nuttall::set(int iSize, bool iPeriodic, float iParam)
{
    mWindow = raisedCosine(
        iSize, iPeriodic, {0.355768, 0.487396, 0.144232, 0.012604}
    );
}

void BlackmanHarris::set(int iSize, bool iPeriodic, float iParam)
{
    mWindow = raisedCosine(
        iSize, iPeriodic, {0.35875, 0.48829, 0.14128, 0.01168}
    );
}

void BlackmanNuttall::set(int iSize, bool iPeriodic, float iParam)
{
    mWindow = raisedCosine(
        iSize, iPeriodic, {0.3635819, 0.4891775, 0.1365995, 0.0106411}
    );
}

void Gaussian::set(int iSize, bool iPeriodic, float iParam)
{
    int d = iPeriodic ? iSize : iSize-1;
    mWindow = var(iSize, 0.0f);
    for (int i=0; i<iSize; i++)
        mWindow[i] = exp(-0.5 * pow( (i-d/2.0) / (iParam*d/2.0), 2 ));
}
