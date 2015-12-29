/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, December 2015
 */

#ifndef COCHLEA_H
#define COCHLEA_H

#include <lube.h>

namespace ssp
{
    class Holdsworth;

    /**
     * Model of a human cochlea; in particular the concept of a filterbank.
     */
    class Cochlea
    {
    public:
        Cochlea(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        ~Cochlea();
        float* operator ()(float iSample, float* oFilter);
    private:
        int mNFilters;
        Holdsworth* mFilter;
        float hzToERB(float iHz);
        float erbToHz(float iERB);
        float hzToERBRate(float iHz);
        float erbRateToHz(float iRate);
        void filter(int iCentre, float iSample, float* oFilter);
    };

    /**
     * One instance of the bandpass filter described by Holdsworth et al. being
     * a component of a gammatone filter bank.
     */
    class Holdsworth
    {
    public:
        Holdsworth();
        void set(float iHz, float iBW, float iPeriod);
        float operator ()(float iSample);
    private:
        static const int cOrder = 4;
        float bwScale();
        float mCentre;
        float mCoeff;
        lube::cfloat mShift;
        lube::cfloat mState[cOrder];
        lube::cfloat mDShift;
        lube::cfloat mUShift;
    };
}

#endif // COCHLEA_H
