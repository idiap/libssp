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

using namespace std;
using namespace ssp;

int main(int argc, char** argv)
{
    // Args
    var arg(argc, argv);

    // Read the waveform from file
    PCM pcm;
    ind fi = arg.index("-f");
    var wav = fi ? arg[fi+1] : "arctic_a0001.wav";
    var a = pcm.read(wav);

    // Find the amax so the y axis can be made symmetric
    var im = lube::iamax(a);
    var m = a[im[0]];

    // Format the plotting commands into a string
    varstream plotstream;
    plotstream
        << "set xrange [0:" << a.size() << "]\n"
        << "set xlabel \"Sample\"\n"
        << "set yrange [" << m * -1 << ":" << m << "]\n"
        << "set ylabel \"Level\"";

    // Plot using gnuplot
    var gnu;
    gnu.push("set grid");
    gnu.push("set style data lines");
    gnu.push(plotstream.str());
    gnu.push("plot \"-\"");
    gnu.push(a);
    file gnuf("gnuplot");
    gnuf.write(lube::nil, gnu);

    // Done
    return 0;
}
