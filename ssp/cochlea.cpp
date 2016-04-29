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


Cochlea::Cochlea()
{
    mNFilters = 0;
}


void Cochlea::set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    // Calculate the centre frequencies equally spaced on ERB rate scale.
    // minHz and maxHz correspond to the range extremes.  We ignore the centres
    // at those frequencies as would a mel filterbank implemented on a DFT.
    mNFilters = iNFilters;
    const bool mel = 0;
    float minRate = mel ? hzToMel(iMinHz) : hzToERBRate(iMinHz);
    float maxRate = mel ? hzToMel(iMaxHz) : hzToERBRate(iMaxHz);
    float step = (maxRate-minRate)/(mNFilters+1);  // +1 = mNFilters+2 filters

    // Go backwards so the cascades work
    for (int i=mNFilters-1; i>=0; --i)
    {
        float rate = minRate + step*(i+1);
        float hz = mel ? melToHz(rate) : erbRateToHz(rate);
        float erb = hzToERB(hz);
        set(i, hz, erb, iPeriod);
    }

    // Zero all the states
    reset();
}

// You'd think there'd be a library function for this
static int factorial(int iN)
{
    int ret = 1;
    while (iN > 1)
        ret *= iN--;
    return ret;
}

float Cochlea::bwScale(int iOrder)
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
    mFilter = 0;
}

Holdsworth::Holdsworth(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    mFilter = new filter[iNFilters];
    Cochlea::set(iMinHz, iMaxHz, iNFilters, iPeriod);
}

Holdsworth::~Holdsworth()
{
    if (mFilter)
        delete [] mFilter;
    mFilter = 0;
}

void Holdsworth::set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    if (mFilter)
        delete [] mFilter;
    mFilter = new filter[iNFilters];
    Cochlea::set(iMinHz, iMaxHz, iNFilters, iPeriod);
}

void Holdsworth::set(int iFilter, float iHz, float iBW, float iPeriod)
{
    filter& f = mFilter[iFilter];
    f.centre = iHz;
    f.coeff = 1.0f - std::exp(-2.0f*PI*iBW/bwScale(cOrder)*iPeriod);
    f.dDelta = std::exp(cfloat(0.0f, -2.0f*PI*iHz*iPeriod));
    f.uDelta = std::exp(cfloat(0.0f,  2.0f*PI*iHz*iPeriod));
    f.dShift = 1.0f;
    f.uShift = 1.0f;
}

void Holdsworth::reset()
{
    for (int i=0; i<mNFilters; i++)
        for (int j=0; j<cOrder+1; j++)
            mFilter[i].state[j] = 0.0f;
}

void Holdsworth::dump()
{
    for (int i=0; i<mNFilters; i++)
    {
        float c = mFilter[i].centre;
        cout << "Centre " << i << ": " << c
             << ", erb: " << hzToERB(c) << endl;
    }
}

void Holdsworth::operator ()(float iSample, float* oFilter)
{
    for (int i=0; i<mNFilters; i++)
    {
        filter& f = mFilter[i];

        // Shift the sample down to DC
        cfloat z = f.dShift * iSample;

        // Filter
        cfloat w;
        for (int j=1; j<cOrder+1; j++)
        {
            w = f.state[j] + f.coeff * (f.state[j-1] - f.state[j]);
            f.state[j-1] = z;
            z = w;
        }
        f.state[cOrder] = z;

        // Shift back up to center frequency
        oFilter[i] = (f.uShift * w).real();

        // Update the frequency shifts
        f.dShift *= f.dDelta;
        f.uShift *= f.uDelta;
    }
}

Lyon::Lyon()
{
    mFilter = 0;
}

Lyon::Lyon(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    mFilter = new filter[iNFilters];
    Cochlea::set(iMinHz, iMaxHz, iNFilters, iPeriod);
}

Lyon::~Lyon()
{
    if (mFilter)
        delete [] mFilter;
    mFilter = 0;
}

void Lyon::set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    if (mFilter)
        delete [] mFilter;
    mFilter = new filter[iNFilters];
    Cochlea::set(iMinHz, iMaxHz, iNFilters, iPeriod);
}

void Lyon::set(int iFilter, float iHz, float iBW, float iPeriod)
{
    filter& f = mFilter[iFilter];
    f.centre = iHz;
    float E = std::exp(-2.0f*PI*iBW/bwScale(cOrder)*iPeriod);
    float C = std::cos( 2.0f*PI*iHz*iPeriod);
    f.coeff[0] = 1.0f - E*C*2 + E*E;
    f.coeff[1] = E*C*2;
    f.coeff[2] = -E*E;
}

void Lyon::reset()
{
    for (int i=0; i<mNFilters; i++)
        for (int j=0; j<cOrder+1; j++)
        {
            mFilter[i].state[0][j] = 0.0f;
            mFilter[i].state[1][j] = 0.0f;
        }
}

void Lyon::dump()
{
    for (int i=0; i<mNFilters; i++)
    {
        float c = mFilter[i].centre;
        cout << "Centre " << i << ": " << c
             << ", erb: " << hzToERB(c) << endl;
    }
}

void Lyon::operator ()(float iSample, float* oFilter)
{
    for (int i=0; i<mNFilters; i++)
    {
        float z = iSample;
        float w;
        filter& f = mFilter[i];
        for (int i=1; i<cOrder+1; i++)
        {
            w = f.coeff[0] * f.state[0][i-1] +
                f.coeff[1] * f.state[0][i] +
                f.coeff[2] * f.state[1][i];
            f.state[1][i-1] = f.state[0][i-1];
            f.state[0][i-1] = z;
            z = w;
        }
        f.state[1][cOrder] = f.state[0][cOrder];
        f.state[0][cOrder] = z;
        oFilter[i] = w;
    }
}


Cascade::Cascade()
{
    mFilter = 0;
}

Cascade::Cascade(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    mFilter = new filter[iNFilters];
    Cochlea::set(iMinHz, iMaxHz, iNFilters, iPeriod);
}

Cascade::~Cascade()
{
    if (mFilter)
        delete [] mFilter;
    mFilter = 0;
}

void Cascade::set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod)
{
    if (mFilter)
        delete [] mFilter;
    mFilter = new filter[iNFilters];
    Cochlea::set(iMinHz, iMaxHz, iNFilters, iPeriod);
}

void Cascade::set(int iFilter, float iHz, float iBW, float iPeriod)
{
    filter& f = mFilter[iFilter];
    f.centre = iHz;
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
    f.filter.set(3, numer, 3, denom);
}

void Cascade::reset()
{
    for (int i=0; i<mNFilters; i++)
        for (int j=0; j<3; j++)
            mFilter[i].state[j] = 0.0f;
}

void Cascade::dump()
{
    for (int i=0; i<mNFilters; i++)
    {
        float c = mFilter[i].centre;
        cout << "Centre " << i << ": " << c
             << ", erb: " << hzToERB(c) << endl;
    }
}

void Cascade::operator ()(float iSample, float* oFilter)
{
    oFilter[mNFilters-1] =
        mFilter[mNFilters-1].filter(iSample, mFilter[mNFilters-1].state);
    for (int i=mNFilters-2; i>=0; --i)
        oFilter[i] = mFilter[i].filter(oFilter[i+1], mFilter[i].state);
}
