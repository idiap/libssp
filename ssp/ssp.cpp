/*
 * Copyright 2014 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, February 2014
 */


#include "ssp.h"

#include <cmath>
#include <random>

using namespace std;
using namespace ssp;

/**
 * Constructor setting sample rate to default (undefined) value.
 */
PCM::PCM(lube::Config& iConfig, var iStr)
    : Config(iConfig)
{
    mStr = iStr;
    mAttr = config();
    mAttr["frameSize"] = 256;
    mAttr["rate"] = 16000;
}


/**
 * Read an audio file
 */
var PCM::read(var iFileName)
{
    var rate = mAttr["rate"].copy();
    file vf("snd", mAttr);
    var snd = vf.read(iFileName);
    if (rate && (rate != mAttr["rate"]))
        throw lube::error("read: file rate does not match pcm rate");
    return snd;
}


/**
 * Write an audio file
 */
void PCM::write(var iFileName, var iVar)
{
    file vf("snd", mAttr);
    vf.write(iFileName, iVar);
}


var PCM::frame(var iVar, int iSize, int iPeriod, bool iPad)
{
    if (iPad)
    {
        // This ensures that frames are aligned in the centre.  Awfully.
        var x(iSize/2, iVar[0]);
        var y(iSize/2, iVar.top());
        x.append(iVar);
        x.append(y);
        iVar = x;
    }

    int nFrames = (iVar.size() - (iSize-iPeriod)) / iPeriod;
    var frame = 0.0f;
    frame.resize(nFrames * iSize);
    var f = frame.view({iSize});
    var g = iVar.view({iSize});
    for (int r=0; r<nFrames; r++)
    {
        f = g;
        if (r != nFrames-1)
        {
            f.advance(iSize);
            g.advance(iPeriod);
        }
    }
    return frame.view({nFrames, iSize});
}

float PCM::hzToRadians(var iHz)
{
    return iHz.cast<float>() / rate() * PI * 2;
}

int PCM::hzToDFTBin(var iHz)
{
    var r = iHz / rate() * mAttr["frameSize"].cast<float>() + 0.5;
    return r.cast<int>();
}

float PCM::dftBinToHz(var iBin)
{
    return rate() * iBin.cast<int>() / mAttr["frameSize"].cast<float>();
}

/**
 * Returns a number of samples corresponding to the given time.  If iPower is
 * AT_LEAST then it is rounded up to the next power of 2; AT_MOST rounds down
 * to the previous power of 2.
 */
int PCM::secondsToSamples(var iSeconds, ind iPower)
{
    int samples = iSeconds.cast<float>() * rate();
    if (!iPower)
        return samples;

    int s = 1;
    while (s < samples)
        s <<= 1;
    if ((s == samples) || (iPower == AT_LEAST))
        return s;
    return s/2;
}

float PCM::samplesToSeconds(var iSamples)
{
    var r = iSamples.cast<float>() / rate();
    return r.cast<float>();
}

Codec::Codec(PCM* iPCM)
{
    mPCM = iPCM;
}

var Codec::encode(var iSignal)
{
    throw lube::error("No encoder defined");
}

var Codec::decode(var iParams)
{
    throw lube::error("No decoder defined");
}

Frame::Frame(int iSize, int iPeriod, bool iPad)
    : UnaryFunctor(iSize)
{
    mPeriod = iPeriod;
    mPad = iPad;
}

var Frame::alloc(var iVar) const
{
    var sh = iVar.shape();
    int n = sh.top().get<int>();
    if (mPad)
        n += mSize;
    int nFrames = (n - (mSize-mPeriod)) / mPeriod;
    sh[sh.size()-1] = nFrames;
    sh.push(mSize);
    var ret = lube::view(sh, iVar.at(0));
    return ret;
}

void Frame::vector(var iVar, var& oVar) const
{
    if (mPad)
    {
        // This ensures that frames are aligned in the centre.  Awfully.
        var x(mSize/2, iVar[0]);
        var y(mSize/2, iVar.top());
        x.append(iVar);
        x.append(y);
        iVar = x;
    }

    var f = oVar.view({mSize});
    var g = iVar.view({mSize});
    int nFrames = oVar.shape(oVar.dim()-2);
    for (int r=0; r<nFrames; r++)
    {
        f = g;
        if (r != nFrames-1)
        {
            f.advance(mSize);
            g.advance(mPeriod);
        }
    }
}


/**
 * General raised cosine window
 */
static var raisedCosine(int iSize, const initializer_list<float> iCoeff)
{
    static const float pi = atan(1.0) * 4;

    var w(iSize, 0.0f);
    for (int i=0; i<iSize; i++)
    {
        int m = -1;
        int j = 0;
        for (const float* c=begin(iCoeff); c!=end(iCoeff); ++c)
        {
            m *= -1;
            w[i] += *c * cos(2.0 * pi * i * j++ / (iSize-1)) * m;
        }
    }
    return w;
}

/**
 * Some particular raised cosine windows
 * http://en.wikipedia.org/wiki/Window_function
 */
var ssp::hanning(int iSize)
{
    return raisedCosine(iSize, {0.5, 0.5});
}

var ssp::hamming(int iSize)
{
    return raisedCosine(iSize, {0.54, 0.46});
}

var ssp::nuttall(int iSize)
{
    return raisedCosine(iSize, {0.355768, 0.487396, 0.144232, 0.012604});
}

var ssp::blackmanharris(int iSize)
{
    return raisedCosine(iSize, {0.35875, 0.48829, 0.14128, 0.01168});
}

var ssp::blackmannuttall(int iSize)
{
    return raisedCosine(iSize, {0.3635819, 0.4891775, 0.1365995, 0.0106411});
}

var ssp::gaussian(int iSize, double iSigma)
{
    var w(iSize, 0.0f);
    for (int i=0; i<iSize; i++)
        w[i] = exp(-0.5 * pow( (i-(iSize-1)/2.0) / (iSigma*(iSize-1)/2.0), 2 ));
    return w;
}

ssp::UnaryFunctor::UnaryFunctor(int iSize)
{
    mDim = 1;
    mSize = iSize;
}

var ssp::UnaryFunctor::alloc(var iVar) const
{
    // Allocate an output array
    var s = iVar.shape();
    s[s.size()-1] = mSize;
    var r = view(s, iVar.at(0));
    return r;
}

ssp::BinaryFunctor::BinaryFunctor(int iSize)
{
    mDim = 1;
    mSize = iSize;
}

var ssp::BinaryFunctor::alloc(var iVar1, var iVar2) const
{
    // Allocate an output array.  Avoid having a trailing dimension of 1 if the
    // size is 1.
    var s = iVar1.shape();
    s[s.size()-1] = mSize;
    if ((s.size() > 1) && (mSize == 1))
        s.pop();
    var r = view(s, iVar1.at(0));
    return r;
}

Autocorrelation::Autocorrelation(int iSize)
    : UnaryFunctor(iSize), mDFT(iSize), mIDFT(iSize, true)
{}

void Autocorrelation::scalar(const var& iVar, var& oVar) const
{
    var tmp = mDFT(iVar);
    tmp.norm();
    tmp /= mSize;
    mIDFT(tmp, oVar);
}

var OverlapAdd::alloc(var iVar) const
{
    if (iVar.dim() < 2)
        throw lube::error("OverlapAdd::alloc: dimension less than 2");
    var sh = iVar.shape();
    int step = sh.pop().get<int>() / 2;
    int size = (sh.top().get<int>() + 1) * step;
    sh[sh.size()-1] = size;
    var ret = lube::view(sh, iVar.at(0));  // Doesn't initialise!  Should it?
    return ret;
}

void OverlapAdd::vector(var iVar, var& oVar) const
{
    // This can surely be vectorised
    int step = iVar.shape(1) / 2;
    for (int i=0; i<oVar.size(); i++)
        oVar(i) = 0.0f;
    for (int f=0; f<iVar.shape(0); f++)
        for (int s=0; s<iVar.shape(1); s++)
            oVar(f*step+s) += iVar(f,s);
}


Filter::Filter(var iNumer, var iDenom)
{
    mDim = 1;
    mNumer = iNumer.copy();
    for(int i=1; i<iDenom.size(); i++)
        mDenom.push(iDenom[i]*-1);
}

void Filter::vector(var iVar, var& oVar) const
{
    // Direct form II filter
    int ns = std::max(mNumer.size(), mDenom.size());
    var state(ns, 0.0f);
    for (int i=0; i<iVar.size(); i++)
    {
        var y = iVar(i);
        if (mDenom)
            y += lube::dot(state, mDenom);
        state.pop();
        state.unshift(y);
        oVar(i) = mNumer ? lube::dot(state, mNumer) : y;
    }
}


var ssp::normal(int iSize, float iMean, float iStdDev)
{
    static std::default_random_engine generator;
    std::normal_distribution<float> distribution(iMean, iStdDev);
    var r(iSize, 0.0f);
    for (int i=0; i<iSize; i++)
        r[i] = distribution(generator);
    return r;
}
