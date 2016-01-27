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
    /** Type of bandpass filter */
    enum {
        BPF_HOLDSWORTH,
        BPF_LYON
    };

    /**
     * Interface to a bandpass filter implementation of a cochlear filter.
     */
    class CochlearBPF
    {
    public:
        virtual void set(float iHz, float iBW, float iPeriod) = 0;
        virtual void reset() = 0;
        virtual float operator ()(float iSample) = 0;
        float centre() { return mCentre; };
    protected:
        float bwScale(int iOrder);
        float mCentre;
    };

    /**
     * Model of a human cochlea; in particular the concept of a filterbank.
     */
    class Cochlea
    {
    public:
        Cochlea(
            float iMinHz, float iMaxHz, int iNFilters, float iPeriod,
            int iType = BPF_HOLDSWORTH
        );
        ~Cochlea();
        float* operator ()(float iSample, float* oFilter);
        void reset();
        void dump();
    private:
        int mNFilters;
        int mType;
        CochlearBPF* mFilter;
        void filter(int iCentre, float iSample, float* oFilter);
    };

    /**
     * One instance of the bandpass filter described by Holdsworth et al. being
     * a component of a gammatone filter bank.
     */
    class Holdsworth : public CochlearBPF
    {
    public:
        Holdsworth();
        void set(float iHz, float iBW, float iPeriod);
        void reset();
        float operator ()(float iSample);
    private:
        static const int cOrder = 4;
        float mCoeff;
        lube::cfloat mShift;
        lube::cfloat mState[cOrder+1];
        lube::cfloat mDShift;
        lube::cfloat mUShift;
    };

    /**
     * One instance of the filter described by Lyon, being a component of an
     * all-pole (non-)gamma-tone filterbank.
     */
    class Lyon : public CochlearBPF
    {
    public:
        Lyon();
        void set(float iHz, float iBW, float iPeriod);
        void reset();
        float operator ()(float iSample);
    private:
        static const int cOrder = 2;
        float mCoeff[3];
        float mState[2][cOrder+1];
    };
}

#endif // COCHLEA_H
