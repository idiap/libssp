/*
 * Copyright 2014 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, July 2014
 */

#ifndef AR_H
#define AR_H

#include <lube.h>
#include "ssp.h"

namespace ssp
{
    int arorder(var iRate);

    /**
     * Functor implementing Levinson-Durbin recursion to calculate
     * autoregression coefficients from autocorrelation.
     */
    class Levinson : public ssp::UnaryFunctor
    {
    public:
        Levinson(int iOrder=0, var iPrior=0.0f);
    protected:
        void vector(
            var iVar, ind iOffsetI, var& oVar, ind iOffsetO
        ) const;
    private:
        var mPrior;
    };

    /**
     * Functor implementing usual AR gain calculation
     */
    class Gain : public ssp::BinaryFunctor
    {
    public:
        Gain(int iOrder=0);
    protected:
        void vector(
            var iVar1, ind iOffset1,
            var iVar2, ind iOffset2,
            var& oVar, ind iOffsetO
        ) const;
    private:
        int mOrder;
    };

    /**
     * Functor to convert AR coefficients to AR spectrum
     */
    class Spectrum : public ssp::BinaryFunctor
    {
    public:
        Spectrum(int iOrder, int iSize);
    protected:
        void vector(
            var iVar1, ind iOffset1,
            var iVar2, ind iOffset2,
            var& oVar, ind iOffsetO
        ) const;
    private:
        int mOrder;
        var mTwiddle;
    };

    /**
     * Retrieve the excitation
     */
    class Excitation : public lube::NaryFunctor
    {
    public:
        Excitation() { mDim = 1; };
    private:
        void vector(var iVar, var& oVar) const;
    };

    /**
     * Resynthesis
     */
    class Resynthesis : public lube::NaryFunctor
    {
    public:
        Resynthesis() { mDim = 1; };
    private:
        void vector(var iVar, var& oVar) const;
    };

    /**
     * Convert AR polynomial to LSPs
     */
    class ToLSP : public ssp::UnaryFunctor
    {
    public:
        ToLSP(int iOrder=0) : ssp::UnaryFunctor(iOrder+2) {};
    private:
        void vector(var iVar, var& oVar) const;
    };

    /**
     * Convert LSPs back to AR polynomial
     */
    class FromLSP : public ssp::UnaryFunctor
    {
    public:
        FromLSP(int iOrder=0) : ssp::UnaryFunctor(iOrder+1) {};
    private:
        void vector(var iVar, var& oVar) const;
    };
}

#endif // AR_H
