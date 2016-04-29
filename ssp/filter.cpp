/*
 * Copyright 2016 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, April 2016
 */

#include <algorithm>
#include <lube/c++blas.h>
#include "filter.h"

using namespace ssp::core;

Filter::Filter()
{
    mNNumer = 0;
    mNumer = 0;
    mNDenom = 0;
    mDenom = 0;
    mNStates = 0;
}

Filter::Filter(int iNNumer, float* iNumer, int iNDenom, float* iDenom)
{
    set(iNNumer, iNumer, iNDenom, iDenom);
}

Filter::~Filter()
{
    if (mNumer)
        delete [] mNumer;
    if (mDenom)
        delete [] mDenom;
    mNumer = 0;
    mDenom = 0;
}

void Filter::set(int iNNumer, float* iNumer, int iNDenom, float* iDenom)
{
    if (iNumer)
    {
        mNNumer = iNNumer;
        mNumer = new float[mNNumer];
        blas::copy(mNNumer, iNumer, mNumer);
    }
    if (iDenom)
    {
        mNDenom = iNDenom-1;
        mDenom = new float[mNDenom];
        blas::copy(mNDenom, iDenom+1, mDenom);
        blas::scal(mNDenom, -1.0f, mDenom);
    }
    mNStates = std::max(mNNumer, mNDenom);
}


/**
 * Filter the input to output.  Allocate the state on the stack.
 */
void Filter::operator()(int iNSamples, float* iSample, float* oSample) const
{
    float state[mNStates];
    for (int i=0; i<mNStates; i++)
        state[i] = 0.0f;
    for (int i=0; i<iNSamples; i++)
        oSample[i] = this->operator()(iSample[i], state);
}

/**
 * Filter one sample at once
 */
float Filter::operator()(float iSample, float* iState) const
{
    float y = iSample;
    if (mDenom)
        y += blas::dot(mNDenom, iState, mDenom);
    for (int s=mNStates-1; s>0; --s)
        iState[s] = iState[s-1];
    iState[0] = y;
    return mNumer ? blas::dot(mNNumer, iState, mNumer) : y;
}
