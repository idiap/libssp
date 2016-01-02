/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, February 2015
 */

#include <cassert>

#include <lube/c++blas.h>
#include "pitch.h"

using namespace ssp;

/**
 * Kalman filter using HNR
 */
class Kalman : public lube::UnaryFunctor
{
public:
    Kalman(var iSeqVar, var iInitMean, var iInitVar);
private:
    void vector(var iVar, var& oVar) const;
    float mSeqVar;
    float mInitMean;
    float mInitVar;
};

Kalman::Kalman(var iSeqVar, var iInitMean, var iInitVar)
{
    mDim = 2;
    mSeqVar = iSeqVar.cast<float>();
    mInitMean = iInitMean.cast<float>();
    mInitVar = iInitVar.cast<float>();
}

#define obs(i) iVar(i,0)
#define obsVar(i) iVar(i,1)
#define stateMean(i) oVar(i,0)
#define stateVar(i) oVar(i,1)
void Kalman::vector(var iVar, var& oVar) const
{
    // Initialise
    stateMean(0) = ( (obs(0) * mInitVar + obsVar(0) * mInitMean) /
                     (obsVar(0) + mInitVar) );
    stateVar(0)  = obsVar(0) * mInitVar / (obsVar(0) + mInitVar);

    // Filter loop
    int n = iVar.shape(0);
    for (int i=1; i<n; i++)
    {
        var predictor = stateVar(i-1) + mSeqVar;
        stateMean(i) = ( (obs(i) * predictor + stateMean(i-1) * obsVar(i)) /
                         (obsVar(i) + predictor) );
        stateVar(i)  = obsVar(i) * predictor / (obsVar(i) + predictor);
    }

    // Smoother loop
    for (int i=n-2; i>=0; --i)
    {
        stateMean(i) = ( stateMean(i+1) * stateVar(i) +
                         stateMean(i)   * mSeqVar );
        stateMean(i) /= (stateVar(i) + mSeqVar);
        var J = stateVar(i) / (stateVar(i) + mSeqVar);
        stateVar(i) = J * (J * stateVar(i+1) + mSeqVar);
    }
}


/**
 * Functor to normalise an autocorrelation by dividing by the first element
 * (the energy).
 */
class NormAC : public lube::UnaryFunctor
{
public:
    NormAC() { mDim = 1; };
private:
    void vector(var iVar, var& oVar) const;
};

void NormAC::vector(var iVar, var& oVar) const
{
    var x = iVar(0);
    x = var(1.0f) / x;
    lube::mul(iVar, x, oVar);
}


/**
 * Finds the index of the maximum value of each array.  Further, ensures that
 * such maxima are not at the borders of the arrays.  It's not all that great;
 * might be better to find the first trough first.
 */
long iamax(int iSize, float* iData, int iLoBin=0, int iHiBin=0)
{
    long lo = iLoBin;
    long hi = iHiBin ? iHiBin : iSize;
    long m = blas::iamax<float>(hi-lo, iData+lo);
    while ( ((m == 0) || (m == hi-lo)) && (lo+1 < hi) )
    {
        if (m == 0)
            lo += 1;
        if (m == hi-lo)
            hi -= 1;
        m = blas::iamax<float>(hi-lo, iData+lo);
    }
    return m+lo;
}


class PitchHNR : public UnaryFunctor
{
public:
    PitchHNR(PCM* iPCM, var iLo, var iHi);
private:
    void vector(var iVar, var& oVar) const;
    PCM* mPCM;
    var mLo;
    var mHi;
    int mLoBin;
    int mHiBin;
};

PitchHNR::PitchHNR(PCM* iPCM, var iLo, var iHi)
    : UnaryFunctor(2)
{
    mLo = iLo;
    mHi = iHi;
    mPCM = iPCM;
    mLoBin = mPCM->secondsToSamples(var(1.0f) / mHi);
    mHiBin = mPCM->secondsToSamples(var(1.0f) / mLo);
};

void PitchHNR::vector(var iVar, var& oVar) const
{
    // The input is nac - normalised autocorrelation
    float* nac = iVar.ptr<float>();
    int size = iVar.size() / 2;
    long pit = iamax(size, nac, mLoBin, mHiBin);
    
    float fnac = std::max(nac[pit], 1e-6f);
    float hnr;
    if ( (nac[pit-1] > nac[pit]) || (nac[pit+1] > nac[pit]) )
        // No peak found; set HNR small
        hnr = 1e-8;
    else
        hnr = fnac / (1.0f - fnac);
    float range = (mHi-mLo).cast<float>();
    oVar[0] = var(1.0f) / mPCM->samplesToSeconds(pit);
    oVar[1] = 1.0f / hnr * range * range;
}

Pitch::Pitch(PCM* iPCM, var iLo, var iHi)
    : UnaryFunctor(2)
{
    mPCM = iPCM;
    mLo = iLo;
    mHi = iHi;
}

void Pitch::scalar(const var& iVar, var& oVar) const
{
    int size = iVar.shape(iVar.dim()-1);
    var w = gaussian(size);
    var f = iVar * w;

    Autocorrelation acorr(size);
    NormAC normac;
    var wac = acorr(w);
    normac(wac, wac);
    for (int i=0; i<wac.size(); i++)
        wac(i) = var(1.0f) / wac(i);

    var ac = acorr(f);
    normac(ac, ac);
    ac *= wac;

    PitchHNR pitchhnr(mPCM, mLo, mHi);
    var prange = mHi-mLo;
    //         SeqVar,       InitMean,       InitVar
    Kalman kalman(1e3, mLo + prange/2, prange*prange);
    var phnr = pitchhnr(ac);
    kalman(phnr, oVar);
#if 0
    // Plot it
    vfile gp("gnuplot");
    var plot;
    int n = oVar.shape(0);
    //plot.push("plot \"-\", \"-\"");
    plot.push("plot \"-\"");
    plot.push(lube::transpose(oVar).view({n}));
    //plot.push(lube::transpose(p).view({n}, n));
    gp.write(lube::nil, plot);
#endif
}

var ssp::excitation(var iPitch, var iHNR, PCM* iPCM)
{
    int framePeriod = 128;
    int frameSize = framePeriod * 2;
    int nFrames = iPitch.size();
    int nSamples = framePeriod * (nFrames-1);

    // Construct an impulse sequence and frame it
    var h(nSamples, 0.0f);
    int i = 0;
    int f = 0;
    while ((i < nSamples) && (f < nFrames))
    {
        int period = iPCM->secondsToSamples(var(1.0) / iPitch[f]);
        if (i + period > nSamples)
            break;
        h[i] = sqrt((float)period);
        i += period;
        f = i / framePeriod;
    }
    Frame frame(frameSize, framePeriod);
    var fh = frame(h);

    // Construct a noise sequence and frame it
    var zero = {-0.5f, 1.0f};
    Filter ri(zero);
    var n = ri(normal(nSamples));
    var fn = frame(n);

    // Add them (into fn)
    assert(fh.size() == fn.size());
    var vn = fn.view({1, frameSize});
    var vh = fh.view({1, frameSize});
    for (f=0; f<nFrames; f++)
    {
        var sn = var(1.0f) / (iHNR[f] + 1.0f);
        var sh = var(1.0f) - sn;
        vn *= sn.sqrt();
        vh *= sh.sqrt();
        if (f < nFrames-1)
        {
            vn.advance(frameSize);
            vh.advance(frameSize);
        }
    }
    fn += fh;

    // Window
    var w = hanning(frameSize+1);
    w.pop();
    fn *= w;

    // Done
    return fn;
}
