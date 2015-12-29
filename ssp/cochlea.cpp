/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, December 2015
 */

#include <cmath>
#include <iostream>
#include "ssp.h"
#include "cochlea.h"

using namespace std;
using namespace ssp;
using namespace lube;


Cochlea::Cochlea(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    // Allocate the filters
    mNFilters = iNFilters;
    mFilter = new Holdsworth[mNFilters];

    // Calculate the centre frequencies equally spaced on ERB rate scale
    float minRate = hzToERBRate(iMinHz);
    float maxRate = hzToERBRate(iMaxHz);
    float step = (maxRate-minRate)/(mNFilters-1);
    for (int i=0; i<mNFilters; i++)
    {
        float rate = minRate + step*i;
        float hz = erbRateToHz(rate);
        float erb = hzToERB(hz);
        cout << "Centre " << i
             << ": " << hz
             << ", erb: " << erb
             << endl;
        mFilter[i].set(hz, erb, iPeriod);
    }
}

Cochlea::~Cochlea()
{
    if (mFilter)
        delete [] mFilter;
    mFilter = 0;
}

float Cochlea::hzToERB(float iHz)
{
    float erb = 24.7f * (4.37e-3*iHz + 1.0f);
    return erb;
}

float Cochlea::erbToHz(float iERB)
{
    float hz = (iERB/24.7f-1.0f) / 4.37e-3;
    return hz;
}

float Cochlea::hzToERBRate(float iHz)
{
    float rate = (1000.0f/107.94f)*std::log(437.0e-3f*iHz+1);
    return rate;
}

float Cochlea::erbRateToHz(float iRate)
{
    float hz = (std::exp(107.94e-3*iRate) - 1.0f) / 437e-3f;
    return hz;
}

float* Cochlea::operator ()(float iSample, float* oFilter)
{
    // Call the filter for each bank; could be parallel.
    for (int f=0; f<mNFilters; f++)
        oFilter[f] = mFilter[f](iSample);
    return oFilter;
}


// You'd think there'd be a library function for this
static int factorial(int iN)
{
    int ret = 1;
    while (iN > 1)
        ret *= iN--;
    return ret;
}

Holdsworth::Holdsworth()
{
    set(0.0f, 0.0f, 0.0f);
}

float Holdsworth::bwScale()
{
    static bool init = false;
    static float a;
    if (!init)
    {
        // Calculate the a_n factor; it's static for a given order
        float numer = PI * factorial(2*cOrder-2) * std::pow(2.0f,-(2*cOrder-2));
        float denom = std::pow(factorial(cOrder-1), 2);
        a = numer / denom;
        init = true;
        cout << "A is: " << a << endl;
    }
    return a;
}

void Holdsworth::set(float iHz, float iBW, float iPeriod)
{
    mCentre = iHz;
    mCoeff = 1.0f - std::exp(-2.0f*PI*iBW/bwScale()*iPeriod);
    mShift = std::exp(cfloat(0.0f, 2.0f*PI*iHz*iPeriod));
    mDShift = 1.0f;
    mUShift = 1.0f;
    for (int j=0; j<cOrder; j++)
        mState[j] = 0.0f;
}

float Holdsworth::operator ()(float iSample)
{
    // Shift the sample down to DC
    cfloat z = mDShift * iSample;

    // Filter
    cfloat w;
    for (int i=1; i<cOrder; i++)
    {
        w = mState[i] + mCoeff * (mState[i-1] - mState[i]);
        mState[i-1] = z;
        z = w;
    }
    mState[cOrder-1] = z;

    // Shift back up to center frequency
    float ret = (mUShift * w).real();

    // Update the frequency shifts
    mDShift *= -mShift;
    mUShift *=  mShift;

    // Return the real result
    return ret;
}
