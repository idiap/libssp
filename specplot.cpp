/*
 * Copyright 2014 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, August 2014
 */

#include <lube.h>
#include "ssp/ssp.h"
#include "ssp/ar.h"

using namespace std;
using namespace ssp;

int main(int argc, char** argv)
{
    // Options
    lube::Option opt("specplot: plot sound spectrogram");
    opt("Waveform file is the first argument");
    opt.parse(argc, argv);

    // Read the waveform from file
    PCM pcm;
    var arg = opt.args();
    if (arg.size() < 1)
        opt.usage(0);
    var wav = arg[0];
    var a = pcm.read(wav);

    // Frame it
    int frameSize = 256;
    int framePeriod = 128;
    var f = pcm.frame(a, frameSize, framePeriod);

    // Window the frames
    var w = nuttall(frameSize);
    f *= w;

    // Choose the output type
    ind ti = arg.index("-t");
    var t = ti ? arg[ti+1] : "spec";
    var p;
    if (t == "spec")
    {
        // Fourier transform
        lube::DFT dft(frameSize);
        var s = dft(f);

        // Periodogram
        p = lube::norm(s);
    }
    else if (t == "ar")
    {
        // LP spectrum
        int order = arorder(pcm.rate());
        Autocorrelation ac(frameSize);
        var af = ac(f);
        Levinson ar(order);
        var lf = ar(af);
        Gain gain(order);
        var g = gain(af, lf);
        Spectrum s(order, 129);
        p = s(lf, g);
    }

    // Plot
    var gnu;
    gnu.push("plot \"-\" matrix with image");
    gnu.push(lube::log(p+1e-8));
    lube::filemodule gnum("gnuplot");
    lube::file& gnuf = gnum.create();
    gnuf.write(lube::nil, gnu);

    return 0;
}
