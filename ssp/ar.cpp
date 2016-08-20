/*
 * Copyright 2014 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, July 2014
 */

#include <cassert>
#include "ar.h"

using namespace ssp;

/**
 * Rule-of-thumb order for AR analysis
 */
int ssp::arorder(var iRate)
{
    int o = int(iRate.cast<float>() / 1000) + 2;
    return o;
}


Levinson::Levinson(int iOrder, var iPrior)
    : ssp::UnaryFunctor(iOrder+1)
{
    mPrior = iPrior;
}


template <class T>
void levinson(int iSize, T iPrior, T* iAC, T* oT)
{
    // Size is order+1
    T c[iSize];
    T p[iSize];
    T* curr = c;
    T* prev = p;
    curr[0] = (T)1.0;
    T error = iAC[0] + iPrior;
    for (int i=1; i<iSize; i++)
    {
        // swap current and previous coefficients
        T* tmp = curr;
        curr = prev;
        prev = tmp;

        // Recurse
        T k = iAC[i];
        for (int j=1; j<i; j++)
            k += prev[j] * iAC[i-j];
        curr[i] = - k / error;
        error *= (T)1.0 - curr[i]*curr[i];
        for (int j=1; j<i; j++)
            curr[j] = prev[j] + curr[i] * prev[i-j];
    }

    // Copy.  Could get rid of this by checking that order is even and swapping
    // curr/prev otherwise.
    for (int i=0; i<iSize; i++)
        oT[i] = curr[i];
}


void Levinson::vector(
    var iVar, ind iOffsetI, var& oVar, ind iOffsetO
) const
{
    switch (iVar.atype())
    {
    case lube::TYPE_FLOAT:
        levinson<float>(mSize, mPrior.cast<float>(),
                        iVar.ptr<float>(iOffsetI),
                        oVar.ptr<float>(iOffsetO)
        );
        break;
    case lube::TYPE_DOUBLE:
        levinson<double>(mSize, mPrior.cast<double>(),
                        iVar.ptr<double>(iOffsetI),
                        oVar.ptr<double>(iOffsetO)
        );
        break;
    default:
        throw std::runtime_error("Levinson::vector: unknown type");
    }
}

Gain::Gain(int iOrder)
    : ssp::BinaryFunctor(1)
{
    mOrder = iOrder;
}

template <class T>
T gain(int iSize, T* iAC, T* iAR)
{
    // This is actually a dot product, i.e., could be optimised
    T g = 0.0;
    for (int i=0; i<iSize; i++)
        g += iAC[i] * iAR[i];
    return g;
}

void Gain::vector(
    var iVar1, ind iOffset1,
    var iVar2, ind iOffset2,
    var& oVar, ind iOffsetO
) const
{
    switch (iVar1.atype())
    {
    case lube::TYPE_FLOAT:
        *oVar.ptr<float>(iOffsetO) = gain<float>(
            mOrder+1, iVar1.ptr<float>(iOffset1), iVar2.ptr<float>(iOffset2)
        );
        break;
    case lube::TYPE_DOUBLE:
        *oVar.ptr<double>(iOffsetO) = gain<double>(
            mOrder+1, iVar1.ptr<double>(iOffset1), iVar2.ptr<double>(iOffset2)
        );
        break;
    default:
        throw std::runtime_error("Gain::vector: unknown type");
    }
}

Spectrum::Spectrum(int iOrder, int iSize)
    : ssp::BinaryFunctor(iSize)
{
    // Set state variables
    mOrder = iOrder;

    // Pre-compute the "twiddle" factors; saves a lot of CPU
    static const float pi = atan(1.0) * 4;
    mTwiddle = lube::view({mSize, mOrder+1}, lube::cfloat(0.0f,0.0f));
    for (int i=0; i<mSize; i++)
        for (int j=0; j<mOrder+1; j++)
            mTwiddle(i,j) =
                lube::exp(lube::cfloat(0.0,-1.0) * pi *
                        (float)i * (float)j / (float)mSize);
}

void Spectrum::vector(
    var iAR, ind iOffsetAR,
    var iGain, ind iOffsetGain,
    var& oVar, ind iOffsetO
) const
{
    if (!iAR.atype<float>())
        throw std::runtime_error("Spectrum::vector: floats only!");

    float* a = iAR.ptr<float>(iOffsetAR);
    float* g = iGain.ptr<float>(iOffsetGain);
    float* o = oVar.ptr<float>(iOffsetO);
    for (int i=0; i<mSize; i++)
    {
        lube::cfloat sm = lube::cfloat(0.0f,0.0f);
        for (int j=0; j<mOrder+1; j++)
            sm += mTwiddle(i,j).get<lube::cfloat>() * a[j];
        float tmp = std::abs(sm);
        o[i] = *g / (tmp*tmp);
    }
}

void Excitation::vector(var iVar, var& oVar) const
{
    // iVar[0]: signal
    // iVar[1]: ar coeffs
    // iVar[2]: gain
    Filter f(iVar[1]);
    f(iVar[0] / lube::sqrt(iVar[2]), oVar);
}

void Resynthesis::vector(var iVar, var& oVar) const
{
    // iVar[0]: excitation
    // iVar[1]: ar coeffs
    // iVar[2]: gain
    Filter f(lube::nil, iVar[1]);
    f(iVar[0] * lube::sqrt(iVar[2]), oVar);
}


/**
 * Convert AR polynomial to LSPs
 *
 * The result includes the values 0 (first) and pi (last).  Of course they are
 * redundant, but this is in the spirit of the first 1 in the polynomial also
 * being redundant.
 */
void ToLSP::vector(var iVar, var& oVar) const
{
    int order = mSize-2;
    var p(mSize, 1.0f);
    var q(mSize, 1.0f);
    q[mSize-1] = -1.0f;
    for (int i=0; i<order; i++)
    {
        p[i+1] = iVar(i+1) + iVar(order-i);
        q[i+1] = iVar(i+1) - iVar(order-i);
    }

    // The roots will not be ordered, but conjugate pairs will be together
    var pr = lube::roots(p);
    var qr = lube::roots(q);

    // From wikipedia:
    // An antipalindromic polynomial (i.e., Q) has 1 as a root.
    // An antipalindromic polynomial of even degree has -1 and 1 as roots
    // A palindromic polynomial (i.e., P) of odd degree has -1 as a root.
    int j = 0;
    for (int i=0; i<mSize-1; i++)
    {
        if (lube::imag(qr[i]) >= 0.0f)
            oVar(j++) = qr[i].arg();
        if (j >= mSize)
        {
            // Some sort of error
            std::cout << "j overflow" << std::endl;
            break;
        }
        if (lube::imag(pr[i]) >= 0.0f)
            oVar(j++) = pr[i].arg();
        if (j >= mSize)
        {
            // Some sort of error
            std::cout << "j overflow" << std::endl;
            break;
        }
    }

    // We can't rely on roots() being ordered, but it's helpful to have 0 and
    // pi at the ends.
    lube::sort(oVar, oVar);
}

/**
 * Convert LSPs back to AR polynomials
 *
 * As for ToLSP, assume that the LSPs are redundant in that they contain 0 and
 * pi.
 */
void FromLSP::vector(var iVar, var& oVar) const
{
    assert(iVar.atype() == lube::TYPE_FLOAT);
    int order = mSize-1;
    var s = lube::sin(iVar);
    var c = lube::cos(iVar);
    var pr(order+2, lube::cfloat(0.0f,0.0f));
    var qr(order+2, lube::cfloat(0.0f,0.0f));

    int qi = 0;
    int pi = 0;
    bool bq = false;

    // Start with (the redundant) iVar(0) = 0
    qr[qi++] = 1.0f;

    // The conjugate pairs alternate
    for (int i=1; i<=order; i++)
    {
        float cf = c[i].get<float>();
        float sf = s[i].get<float>();
        if (bq)
        {
            qr[qi++] = lube::cfloat(cf,  sf);
            qr[qi++] = lube::cfloat(cf, -sf);
        }
        else
        {
            pr[pi++] = lube::cfloat(cf,  sf);
            pr[pi++] = lube::cfloat(cf, -sf);
        }
        bq = !bq;
    }

    // The location of the remaining root depends on the order
    if (bq)
        qr[qi] = -1.0f;
    else
        pr[pi] = -1.0f;
    var p = lube::poly(pr);
    var q = lube::poly(qr);
    q += p;
    lube::real(q, oVar);
    oVar *= 0.5f;
}
