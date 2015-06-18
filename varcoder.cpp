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
#include "arcodec.h"

using namespace std;
using namespace ssp;

void usage()
{
    cout << "varcoder" << endl;
}

int main(int argc, char** argv)
{
    // What to do
    bool encode = false;
    bool decode = false;
    bool oracle = false;

    // Command line
    lube::Config cnf;
    lube::Option opt(argc, argv, "deoC:");
    while (opt)
    {
        file ini("ini");
        switch (opt.get())
        {
        case 'd':
            decode = true;
            break;
        case 'e':
            encode = true;
            break;
        case 'o':
            oracle = true;
            break;
        case 'C':
            cnf.read(opt.arg());
            break;
        default:
            usage();
            exit(EXIT_FAILURE);
        }
    }

    // Register rest of the command line
    var arg = opt;
    if (argc < 3)
        throw lube::error("Not enough args");

    // The input and output should be the last two arguments
    var ofile = arg.pop();
    var ifile = arg.pop();

    // AR codec
    PCM pcm(cnf);
    ARCodec arcodec(&pcm, oracle);
    file htk("htk");

    if (!encode && !decode)
    {
        // Best effort copy synthesis
        var a = pcm.read(ifile);
        var params = arcodec.encode(a);
        var signal = arcodec.decode(params);
        pcm.write(ofile, signal);
    }

    if (encode)
    {
        // Read the file, hanging on to the attributes
        var a = pcm.read(ifile);
        var params = arcodec.encode(a);
        arcodec.write(ofile, params);
    }

    if (decode)
    {
        // Read the params
        var params = arcodec.read(ifile);
        var signal = arcodec.decode(params);
        pcm.write(ofile, signal);
    }

    // Done
    return 0;
}
