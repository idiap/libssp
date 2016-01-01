/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, December 2015
 */

#include <cmath>
#include <cassert>
#include <iostream>
#include "ssp.h"
#include "cochlea.h"

using namespace std;
using namespace ssp;
using namespace lube;


Cochlea::Cochlea(
    float iMinHz, float iMaxHz, int iNFilters, float iPeriod, int iType
)
{
    // Allocate the filters
    mNFilters = iNFilters;
    mType = iType;
    switch (mType)
    {
    case BPF_HOLDSWORTH:
        mFilter = new Holdsworth[mNFilters];
        break;
    case BPF_LYON:
        mFilter = new Lyon[mNFilters];
        break;
    default:
        assert(0);
    }

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
        switch (mType)
        {
        case BPF_HOLDSWORTH:
            dynamic_cast<Holdsworth*>(mFilter)[i].set(hz, erb, iPeriod);
            break;
        case BPF_LYON:
            dynamic_cast<Lyon*>(mFilter)[i].set(hz, erb, iPeriod);
            break;
        default:
            assert(0);
        }
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
    switch (mType)
    {
    case BPF_HOLDSWORTH:
        for (int f=0; f<mNFilters; f++)
            oFilter[f] = dynamic_cast<Holdsworth*>(mFilter)[f](iSample);
        break;
    case BPF_LYON:
        for (int f=0; f<mNFilters; f++)
            oFilter[f] = dynamic_cast<Lyon*>(mFilter)[f](iSample);
        break;
    default:
        assert(0);
    }
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

float CochlearBPF::bwScale(int iOrder)
{
    // This will break if two types of filter use it for different orders in
    // the same instance
    static bool init = false;
    static float a;
    if (!init)
    {
        // Calculate the a_n factor; it's static for a given order
        float numer = PI * factorial(2*iOrder-2) * std::pow(2.0f,-(2*iOrder-2));
        float denom = std::pow(factorial(iOrder-1), 2);
        a = numer / denom;
        init = true;
        cout << "A is: " << a << endl;
    }
    return a;
}

Holdsworth::Holdsworth()
{
    set(0.0f, 0.0f, 0.0f);
}

void Holdsworth::set(float iHz, float iBW, float iPeriod)
{
    mCentre = iHz;
    mCoeff = 1.0f - std::exp(-2.0f*PI*iBW/bwScale(cOrder)*iPeriod);
    mShift = std::exp(cfloat(0.0f, 2.0f*PI*iHz*iPeriod));
    mDShift = 1.0f;
    mUShift = 1.0f;
    for (int j=0; j<cOrder+1; j++)
        mState[j] = 0.0f;
}

float Holdsworth::operator ()(float iSample)
{
    // Shift the sample down to DC
    cfloat z = mDShift * iSample;

    // Filter
    cfloat w;
    for (int i=1; i<cOrder+1; i++)
    {
        w = mState[i] + mCoeff * (mState[i-1] - mState[i]);
        mState[i-1] = z;
        z = w;
    }
    mState[cOrder] = z;

    // Shift back up to center frequency
    float ret = (mUShift * w).real();

    // Update the frequency shifts
    mDShift *= -mShift;
    mUShift *=  mShift;

    // Return the real result
    return ret;
}

Lyon::Lyon()
{
    set(0.0f, 0.0f, 0.0f);
}

void Lyon::set(float iHz, float iBW, float iPeriod)
{
    mCentre = iHz;
    float e = std::exp(-2.0f*PI*iBW/bwScale(cOrder)*iPeriod);
    float s = std::sin(2.0f*PI*iHz*iPeriod);
    float c = std::cos(2.0f*PI*iHz*iPeriod);
    mCoeff[0] = e*s;
    mCoeff[1] = e*c*2;
    mCoeff[2] = -e*e;
    for (int j=0; j<cOrder+1; j++)
    {
        mState[0][j] = 0.0f;
        mState[1][j] = 0.0f;
    }
}

float Lyon::operator ()(float iSample)
{
    // Filter
    float z = iSample;
    float w;
    for (int i=1; i<cOrder+1; i++)
    {
        w = mCoeff[0] * mState[0][i-1] +
            mCoeff[1] * mState[0][i] +
            mCoeff[2] * mState[1][i];
        mState[1][i-1] = mState[0][i-1];
        mState[0][i-1] = z;
        z = w;
    }
    mState[1][cOrder] = mState[0][cOrder];
    mState[0][cOrder] = z;

    // Return the real result
    return w;
}
