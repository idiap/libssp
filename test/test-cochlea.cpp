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
    var noise = normal(nSamples, 0.0f, 1.0f);
    var filter = lube::view({nSamples, nFilters}, 0.0f);
    float* n = noise.ptr<float>();
    float* f = filter.ptr<float>();
    for (int s=0; s<nSamples; s++)
        c(n[s], &f[s*nFilters]);

    filter.transpose();
    cout << "Shape: " << filter.shape() << endl;
    lube::DFT dft(nSamples, 0.0f);
    var t = dft(filter);
    var d = lube::norm(t).log() / lube::log(10.0f) * 10.0f;
    cout << "Norm: " << d.shape() << endl;

    // XScale
    int nDFT = nSamples/2+1;
    float step = rate / 2 / (nDFT-1);
    var xscale = lube::range(0.0f, rate/2, step);
    
    // Plot it; ...which turns out to be a bit of a faff
    var gnu;
    gnu.push("set grid");
    gnu.push("set nokey");
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
