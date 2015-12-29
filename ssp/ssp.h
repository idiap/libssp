/*
 * Copyright 2014 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, February 2014
 */

#ifndef SSP_H
#define SSP_H


#include <lube.h>
#include <lube/config.h>


namespace ssp
{
    // Constants
    const double PI = lube::atan(1.0).cast<double>()*4;

    /**
     * PCM class
     *
     * A configurable object that knows about various rates and conversions.
     * Defined by at least a sample rate.
     */
    class PCM : public lube::Config
    {
    public:
        PCM(var iStr="PCM") { mStr = iStr; };
        PCM(lube::Config& iConfig, var iStr="PCM");

        enum {
            AT_LEAST,
            AT_MOST
        };

        var read(var iFileName);
        void write(var iFileName, var iVar);
        var frame(var iVar, int iSize, int iPeriod, bool iPad=true);
        float rate() { return mAttr["rate"].cast<float>(); };
        float hzToRadians(var iHz);
        int hzToDFTBin(var iHz);
        float dftBinToHz(var iBin);
        int secondsToSamples(var iSeconds, ind iPower=-1);
        float samplesToSeconds(var iSamples);
    private:
        var mAttr;
    };


    /**
     * Abstract codec class
     */
    class Codec
    {
    public:
        Codec(PCM* iPCM);
        virtual ~Codec() {};
        virtual var encode(var iSignal);
        virtual var decode(var iParams);
    protected:
        PCM* mPCM;
    };


    /*
     * Windows
     */
    var hanning(int iSize);
    var hamming(int iSize);
    var nuttall(int iSize);
    var blackmanharris(int iSize);
    var blackmannuttall(int iSize);
    var gaussian(int iSize, double iSigma=0.5);

    /*
     * Functors
     */
    class UnaryFunctor : public lube::UnaryFunctor
    {
    public:
        UnaryFunctor(int iSize);
    protected:
        var alloc(var iVar) const;
        int mSize;
    };

    class BinaryFunctor : public lube::BinaryFunctor
    {
    public:
        BinaryFunctor(int iSize);
    protected:
        var alloc(var iVar1, var iVar2) const;
        int mSize;
    };


    class Autocorrelation : public UnaryFunctor
    {
    public:
        Autocorrelation(int iSize);
    protected:
        void scalar(const var& iVar, var& oVar) const;
    private:
        lube::DFT mDFT;
        lube::DFT mIDFT;
    };

    class Frame : public UnaryFunctor
    {
    public:
        Frame(int iSize, int iPeriod, bool iPad=true);
    protected:
        var alloc(var iVar) const;
        void vector(var iVar, var& oVar) const;
    private:
        int mPeriod;
        bool mPad;
    };

    class OverlapAdd : public lube::UnaryFunctor
    {
    public:
        OverlapAdd() { mDim = 2; };
    protected:
        var alloc(var iVar) const;
        void vector(var iVar, var& oVar) const;
    };

    class Filter : public lube::UnaryFunctor
    {
    public:
        Filter(var iNumer, var iDenom=lube::nil);
    protected:
        void vector(var iVar, var& oVar) const;
    private:
        var mNumer;
        var mDenom;
    };

    var normal(int iSize, float iMean=0.0f, float iStdDev=1.0f);
}


#endif // SSP_H
