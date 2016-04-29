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
#include "filter.h"

namespace ssp
{
    /** Type of bandpass filter */
    enum {
        BPF_HOLDSWORTH,
        BPF_LYON,
        BPF_CASCADE
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
        lube::cfloat mDDelta;
        lube::cfloat mUDelta;
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

    /**
     * One instance of the cascade filter described by Lyon, being a component
     * of a two-pole two-zero filterbank.
     */
    class Cascade : public CochlearBPF
    {
    public:
        Cascade();
        void set(float iHz, float iBW, float iPeriod);
        void reset();
        float operator ()(float iSample);
    private:
        core::Filter mFilter;
        float mState[3];
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
        union {
            // Polymorphism ought to work here, but it leads to lots of
            // dynamic_casts because it's an array.
            Holdsworth* holdsworth;
            Lyon* lyon;
            Cascade* cascade;
        } mFilter;
        void filter(int iCentre, float iSample, float* oFilter);
    };
}

#endif // COCHLEA_H
