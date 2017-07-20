/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, February 2015
 */

#include <lube.h>
#include <lube/config.h>
#include "ssp/arcodec.h"

using namespace std;
using namespace ssp;

int main(int argc, char** argv)
{
    // Command line
    lube::Option opt("varcoder: SSP based vocoder");
    opt(" <args> should be the input and output files respectively");
    opt('e', "Read a wave file and encode as parameters");
    opt('d', "Read parameters and decode");
    opt('o', "Use the oracle excitation in the AR codec");
    opt('C', "Read configuration file", "/dev/null");
    opt("Default behaviour is a best-effort encode-decode copy");
    opt.parse(argc, argv);

    // Configuration
    lube::Config cnf;
    if (opt['C'] != "/dev/null")
        cnf.configFile(opt['C']);

    // Register rest of the command line
    var arg = opt.args();
    if (arg.size() < 2)
        opt.usage(0);

    // The input and output should be the last two arguments
    var ofile = arg.pop();
    var ifile = arg.pop();

    // AR codec
    PCM pcm;
    ARCodec arcodec(&pcm, bool(opt['o']));

    if (!opt['e'] && !opt['d'])
    {
        // Best effort copy synthesis
        var a = pcm.read(ifile);
        var params = arcodec.encode(a);
        var signal = arcodec.decode(params);
        pcm.write(ofile, signal);
    }

    if (opt['e'])
    {
        // Read the file, hanging on to the attributes
        var a = pcm.read(ifile);
        var params = arcodec.encode(a);
        arcodec.write(ofile, params);
    }

    if (opt['d'])
    {
        // Read the params
        var params = arcodec.read(ifile);
        var signal = arcodec.decode(params);
        pcm.write(ofile, signal);
    }

    // Done
    return 0;
}
