/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, May 2015
 */

#include "arcodec.h"
#include "ar.h"
#include "pitch.h"

using namespace ssp;

ARCodec::ARCodec(PCM* iPCM, bool iOracle)
    : Codec(iPCM)
{
    mOracle = iOracle;
}

var ARCodec::encode(var iSignal)
{
    // Frame and window.  The window should be asymmetric, so ask for one too
    // long, then pop off the last sample.
    int framePeriod = mPCM->samplesFromSeconds(0.005, PCM::AT_LEAST);
    int frameSize = framePeriod * 2;
    int pitchSize = mPCM->samplesFromSeconds(0.025, PCM::AT_LEAST);
    Frame frame(frameSize, framePeriod);
    var f = frame(iSignal);
    var w = hanning(frameSize+1);
    w.pop();
    f *= w;

    // AR analysis
    int order = arorder(mPCM->rate());
    Autocorrelation acorr(frameSize);
    Levinson lev(order);
    Gain gain(order);
    ToLSP toLSP(order);

    var ac = acorr(f);
    var ar = lev(ac);
    var gg = gain(ac, ar);

    var ret;
    ret[0] = toLSP(ar);
    ret[1] = gg;

    if (mOracle)
    {
        // Oracle excitation
        Excitation excit;
        ret[2] = excit({f, ar, gg});
    }
    else
    {
        // Pitch / HNR excitation
        Frame pframe(pitchSize, framePeriod);
        var pf = pframe(iSignal);
        Pitch pitch(mPCM);
        var p = pitch(pf);

        // pitch & hnr should be separate
        p.transpose();
        ret[2] = p.subview(1, 0);
        ret[3] = p.subview(1, ret[2].size());
    }

    return ret;
}

var ARCodec::decode(var iParams)
{
    // Resynthesise using decomposed signals
    int order = arorder(mPCM->rate());
    Resynthesis resynth;
    FromLSP fromLSP(order);

    var ar = fromLSP(iParams[0]);
    var re;
    if (mOracle)
        re = resynth({iParams[2],ar,iParams[1]});
    else
    {
        var ex = excitation(iParams[2], iParams[3], mPCM);
        re = resynth({ex, ar, iParams[1]});
    }

    // Reconstruction using overlap-add
    OverlapAdd ola;
    var o = ola(re);

    return o;
}

var ARCodec::read(var iFile)
{
    // Begin by reading the HTK file, which is LSPs and log(gg)
    file htk("htk");
    var prmFile = iFile;
    var prm = htk.read(prmFile);

    // Split it into the right parts
    int nFrames = prm.shape(0);
    int nParams = prm.shape(1);
    var lsp = lube::view({nFrames,nParams+1});
    var lgg(nFrames, 0.0f);
    for (int f=0; f<nFrames; f++)
    {
        static const float pi = atan(1.0) * 4;
        lsp(f, 0) = 0.0f;
        for (int p=0; p<nParams-1; p++)
            lsp(f, p+1) = prm(f, p);
        lsp(f, nParams) = pi;
        lgg(f) = prm(f, nParams-1);
    }

    // The log(f0) and hnr parts are text
    file txt("txt");
    var lf0File = iFile.replace("\\....", ".lf0");
    var hnrFile = iFile.replace("\\....", ".hnr");
    var lf0 = txt.read(lf0File);
    var tmp = txt.read(hnrFile);

    // The library should do this as a broadcast cast<>()!
    var f0(nFrames, 0.0f);
    var hnr(nFrames, 0.0f);
    for (int i=0; i<nFrames; i++)
    {
        f0[i] = lube::exp(lf0[i].cast<float>());
        hnr[i] = tmp[i].cast<float>();
    }

    // Return a "tuple"
    return {lsp, lube::exp(lgg), f0, hnr};
}

void ARCodec::write(var iFile, var iParams)
{
    // We need to concatenate the meaningful parts of the LSP and log(gg)
    // together into one file
    var lsp = iParams[0];
    var lgg = lube::log(iParams[1]);
    var lf0 = lube::log(iParams[2]);
    var hnr = iParams[3];

    int nFrames = lsp.shape(0);
    int nParams = lsp.shape(1)-1;
    var prm = lube::view({nFrames, nParams});
    for (int f=0; f<nFrames; f++)
    {
        for (int p=0; p<nParams-1; p++)
            prm(f, p) = lsp(f, p+1);
        prm(f, nParams-1) = lgg[f];
    }

    // Params are HTK format
    file htk("htk");
    var prmFile = iFile;
    htk.write(prmFile, prm);

    // log(f0) and hnr files are text
    file txt("txt");
    var lf0File = iFile.replace("\\....", ".lf0");
    var hnrFile = iFile.replace("\\....", ".hnr");
    txt.write(lf0File, lf0);
    txt.write(hnrFile, hnr);
}
