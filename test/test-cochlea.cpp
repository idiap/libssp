/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, December 2015
 */

#include <lube.h>
#include "ssp/ssp.h"
#include "ssp/cochlea.h"

using namespace std;
using namespace ssp;

int main(int argc, char** argv)
{
    // Set the FP precision to be less than the difference between different
    // numerical libraries
    // std::cout.precision(4);

    // Some parameters
    float rate = 16000;
    float period = 1.0f/rate;
    int nFilters = 10;
    int nSamples = 2048;

    // Stick noise through a cochlea
    Cochlea c(100, rate/2*0.8f, nFilters, period);
    lube::DFT dft(nSamples, 0.0f);
#if 1
    // Use noise; a window will be necessary
    int nAverage = 100;
    var filter = lube::view({nFilters, nSamples/2+1}, 0.0f);
    var w = nuttall(nSamples);
    filter *= 0.0f;
    cout << "Filter: " << filter.shape() << endl;
    for (int i=0; i<nAverage; i++)
    {
        var noise = normal(nSamples, 0.0f, 1.0f) / sqrt(nSamples);
        noise *= w;
        var nfilt = lube::view({nSamples, nFilters}, 0.0f);
        nfilt *= 0.0f;
        float* n = noise.ptr<float>();
        float* f = nfilt.ptr<float>();
        for (int s=0; s<nSamples; s++)
            c(n[s], &f[s*nFilters]);

        // Filter rows are channels; transpose so they're samples and DFT
        nfilt.transpose();
        var t = dft(nfilt);
        var d = lube::norm(t);
        filter += d;
    }
    var d = (filter / nAverage).log() / lube::log(10.0f) * 10.0f;
#else
    // Use impulse; no window is necessary
    var impulse(nSamples, 0.0f);
    impulse[0] = 1.0f;
    var ifilter = lube::view({nSamples, nFilters}, 0.0f);
    ifilter *= 0.0f;
    float* n = impulse.ptr<float>();
    float* f = ifilter.ptr<float>();
    for (int s=0; s<nSamples; s++)
        c(n[s], &f[s*nFilters]);

    // Filter rows are channels; transpose so they're samples and DFT
    ifilter.transpose();
    var t = dft(ifilter);
    var d = lube::norm(t).log() / lube::log(10.0f) * 10.0f;
#endif
    cout << "Norm: " << d.shape() << endl;

    // XScale
    int nDFT = nSamples/2+1;
    float step = rate / 2 / (nDFT-1);
    var xscale = lube::range(0.0f, rate/2, step);
    
    // Plot it; ...which turns out to be a bit of a faff
    var gnu;
    gnu.push("set grid");
    gnu.push("set xlabel \"Frequency\"");
    gnu.push("set ylabel \"dB\"");
    gnu.push("set nokey");
#if 0
    gnu.push("set logscale x");
    gnu.push("set xrange [30:8000]");
#endif
    gnu.push("set style data lines");
    varstream vs;
    vs << "plot \"-\" using 1:2";
    for (int i=1; i<nFilters; i++)
        vs << ", \"-\" using 1:2";
    gnu.push(var(vs));
    for (int f=0; f<nFilters; f++)
    {
        var data = lube::view({2, nDFT}, 0.0f);
        for (int i=0; i<nDFT; i++)
        {
            data(0,i) = xscale[i];
            data(1,i) = d(f,i);
        }
        gnu.push(data);
    }
    file gnuf("gnuplot");
    // gnuf.write(lube::nil, gnu); // For immediate display
    gnuf.write("test-cochlea.eps", gnu); // For file

    // Done
    return 0;
}
