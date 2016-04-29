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
    /**
     * Model of a human cochlea; in particular the concept of a filterbank.
     */
    class Cochlea
    {
    public:
        Cochlea();
        virtual ~Cochlea() {};
        void set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        virtual void operator ()(float iSample, float* oFilter) = 0;
        virtual void reset() = 0;
        virtual void dump() = 0;
    protected:
        virtual void set(int iFilter, float iHz, float iBW, float iPeriod) = 0;
        float bwScale(int iOrder);
        int mNFilters;
    };

    /**
     * The gammatone filterbank described by Holdsworth et al.
     */
    class Holdsworth : public Cochlea
    {
    public:
        Holdsworth();
        Holdsworth(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        ~Holdsworth();
        void set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        void reset();
        void dump();
        void operator ()(float iSample, float* oFilter);
    protected:
        void set(int iFilter, float iHz, float iBW, float iPeriod);
    private:
        static const int cOrder = 4;
        struct filter
        {
            float centre;
            float coeff;
            lube::cfloat dDelta;
            lube::cfloat uDelta;
            lube::cfloat state[cOrder+1];
            lube::cfloat dShift;
            lube::cfloat uShift;
        };
        filter* mFilter;
    };

    /**
     * The all-pole (non-)gamma-tone filterbank of Lyon.
     */
    class Lyon : public Cochlea
    {
    public:
        Lyon();
        ~Lyon();
        Lyon(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        void set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        void reset();
        void dump();
        void operator ()(float iSample, float* oFilter);
    protected:
        void set(int iFilter, float iHz, float iBW, float iPeriod);
    private:
        static const int cOrder = 2;
        struct filter
        {
            float centre;
            float coeff[3];
            float state[2][cOrder+1];
        };
        filter* mFilter;
    };

    /**
     * Lyon's two-pole two-zero cascade filterbank.
     */
    class Cascade : public Cochlea
    {
    public:
        Cascade();
        ~Cascade();
        Cascade(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        void set(float iMinHz, float iMaxHz, int iNFilters, float iPeriod);
        void reset();
        void dump();
        void operator ()(float iSample, float* oFilter);
    protected:
        void set(int iFilter, float iHz, float iBW, float iPeriod);
    private:
        struct filter
        {
            float centre;
            core::Filter filter;
            float state[3];
        };
        filter* mFilter;
    };
}

#endif // COCHLEA_H
