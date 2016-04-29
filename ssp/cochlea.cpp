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
#include "warp.h"

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
        mFilter.holdsworth = new Holdsworth[mNFilters];
        break;
    case BPF_LYON:
        mFilter.lyon = new Lyon[mNFilters];
        break;
    case BPF_CASCADE:
        mFilter.cascade = new Cascade[mNFilters];
        break;
    default:
        assert(0);
    }

    // Calculate the centre frequencies equally spaced on ERB rate scale.
    // minHz and maxHz correspond to the range extremes.  We ignore the centres
    // at those frequencies as would a mel filterbank implemented on a DFT.
    const bool mel = 0;
    float minRate = mel ? hzToMel(iMinHz) : hzToERBRate(iMinHz);
    float maxRate = mel ? hzToMel(iMaxHz) : hzToERBRate(iMaxHz);
    float step = (maxRate-minRate)/(mNFilters+1);  // +1 = mNFilters+2 filters
    for (int i=0; i<mNFilters; i++)
    {
        float rate = minRate + step*(i+1);
        float hz = mel ? melToHz(rate) : erbRateToHz(rate);
        float erb = hzToERB(hz);
        switch (mType)
        {
        case BPF_HOLDSWORTH:
            mFilter.holdsworth[i].set(hz, erb, iPeriod);
            break;
        case BPF_LYON:
            mFilter.lyon[i].set(hz, erb, iPeriod);
            break;
        case BPF_CASCADE:
            mFilter.cascade[i].set(hz, erb, iPeriod);
            break;
        default:
            assert(0);
        }
    }
}

Cochlea::~Cochlea()
{
    // Any type will do
    // ...no, it won't
    if (mFilter.cascade)
        delete [] mFilter.cascade;
    mFilter.lyon = 0;
}

float* Cochlea::operator ()(float iSample, float* oFilter)
{
    // Call the filter for each bank; could be parallel.
    switch (mType)
    {
    case BPF_HOLDSWORTH:
        for (int f=0; f<mNFilters; f++)
            oFilter[f] = mFilter.holdsworth[f](iSample);
        break;
    case BPF_LYON:
        for (int f=0; f<mNFilters; f++)
            oFilter[f] = mFilter.lyon[f](iSample);
        break;
    case BPF_CASCADE:
        oFilter[mNFilters-1] = mFilter.cascade[mNFilters-1](iSample);
        for (int f=mNFilters-2; f>=0; --f)
            oFilter[f] = mFilter.cascade[f](oFilter[f+1]);
        break;
    default:
        assert(0);
    }
    return oFilter;
}


void Cochlea::reset()
{
    switch (mType)
    {
    case BPF_HOLDSWORTH:
        for (int f=0; f<mNFilters; f++)
            mFilter.holdsworth[f].reset();
        break;
    case BPF_LYON:
        for (int f=0; f<mNFilters; f++)
            mFilter.lyon[f].reset();
        break;
    case BPF_CASCADE:
        for (int f=0; f<mNFilters; f++)
            mFilter.cascade[f].reset();
        break;
    default:
        assert(0);
    }
}

void Cochlea::dump()
{
    for (int i=0; i<mNFilters; i++)
    {
        cout << "Centre " << i << ": ";
        float c;
        switch (mType)
        {
        case BPF_HOLDSWORTH:
            c = mFilter.holdsworth[i].centre();
            break;
        case BPF_LYON:
            c = mFilter.lyon[i].centre();
            break;
        case BPF_CASCADE:
            c = mFilter.cascade[i].centre();
            break;
        default:
            assert(0);
        }
        cout << c << ", erb: " << hzToERB(c) << endl;
    }
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
    mDDelta = std::exp(cfloat(0.0f, -2.0f*PI*iHz*iPeriod));
    mUDelta = std::exp(cfloat(0.0f,  2.0f*PI*iHz*iPeriod));
    mDShift = 1.0f;
    mUShift = 1.0f;
    reset();
}

void Holdsworth::reset()
{
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
    mDShift *= mDDelta;
    mUShift *= mUDelta;

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
    float c = std::cos(2.0f*PI*iHz*iPeriod);
    mCoeff[0] = 1.0f - e*c*2 + e*e;
    mCoeff[1] = e*c*2;
    mCoeff[2] = -e*e;
    reset();
}

void Lyon::reset()
{
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

    // Return the result
    return w;
}


Cascade::Cascade()
{
    set(0.0f, 0.0f, 0.0f);
}

void Cascade::set(float iHz, float iBW, float iPeriod)
{
    mCentre = iHz;
    float Ep = std::exp(-2.0f*PI*iBW/bwScale(1)*iPeriod);
    float Cp = std::cos(2.0f*PI*iHz*iPeriod);
    float Ez = std::exp(-2.0f*PI*iBW/bwScale(1)*2.0*iPeriod);
    float Cz = std::cos(2.0f*PI*iHz*iPeriod*1.5);
    float numer[3];
    float denom[3];
    float A = ( (1.0f - Ep*Cp*2 + Ep*Ep) /
                (1.0f - Ez*Cz*2 + Ez*Ez) );
    numer[0] =  A;
    numer[1] = -A*Ez*Cz*2;
    numer[2] =  A*Ez*Ez;
    denom[0] =  1.0f;
    denom[1] = -Ep*Cp*2;
    denom[2] =  Ep*Ep;
    mFilter.set(3, numer, 3, denom);
    reset();
}

void Cascade::reset()
{
    for (int j=0; j<3; j++)
        mState[j] = 0.0f;
}

float Cascade::operator ()(float iSample)
{
    // Filter
    float w = mFilter(iSample, mState);

    // Return the result
    return w;
}
