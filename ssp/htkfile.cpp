/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, April 2015
 */

#include <cassert>
#include <fstream>
#include <lube/module.h>

namespace libube
{
    /**
     * HTK parameter kinds
     */
    enum {
        HTK_LPC =       1,
        HTK_LPREFC =    2,
        HTK_LPCEPSTRA = 3,
        HTK_MFCC =      6,
        HTK_FBANK =     7,
        HTK_MELSPEC =   8,
        HTK_USER =      9,
        HTK_PLP =      11,
        HTK_E =   0000100,
        HTK_N =   0000200,
        HTK_D =   0000400,
        HTK_A =   0001000,
        HTK_Z =   0004000,
        HTK_0 =   0020000,
        HTK_T =   0100000
    };

    /**
     * Class to handle HTK feature format files
     */
    class HTK : public file
    {
    public:
        HTK(var iAttr);
        virtual var read(var iFile);
        virtual void write(var iFile, var iVar);

    private:
        var mAttr;
    };


    void factory(Module** oModule, var iArg)
    {
        *oModule = new HTK(iArg);
    }
}


using namespace libube;

HTK::HTK(var iAttr)
{
    // Ring alarm bells
    assert(sizeof(float) == 4);
    assert(sizeof(int)   == 4);
    assert(sizeof(short) == 2);

    // Set up the attributes
    mAttr["kind"] = HTK_USER;
    mAttr["period"] = 0.01f;
    for (int i=0; i<iAttr.size(); i++)
        mAttr[iAttr.key(i)] = iAttr[i];
};

var HTK::read(var iFile)
{
    // Open the file as a stream
    std::ifstream is(iFile.str(), std::ifstream::in | std::ifstream::binary);
    if (is.fail())
        throw error("htkfile::read(): Open failed");

    // The 12 byte header
    int nSamples;
    int sampPeriod;
    short sampSize;
    short parmKind;
    
    // Read the 12 byte header
    is.read((char*)&nSamples, 4);
    is.read((char*)&sampPeriod, 4);
    is.read((char*)&sampSize, 2);
    is.read((char*)&parmKind, 2);

    // Store metadata
    mAttr["kind"] = parmKind;
    mAttr["period"] = (float)sampPeriod * 1e-7;

    // And the rest
    int n = nSamples * sampSize;
    var data(n, 0.0f);
    is.read(data.ptr<char>(), n);
    int size = sampSize / sizeof(float);
    return data.view({nSamples, size});
}

void HTK::write(var iFile, var iVar)
{
    // Prepare the header
    int nSamples = iVar.shape(-2);
    int sampPeriod = mAttr["period"].cast<float>() * 1e7 + 0.5;
    short sampSize = iVar.shape(-1) * sizeof(float);
    short parmKind = mAttr["kind"].cast<int>();

    // Open the file
    std::ofstream os(iFile.str(), std::ofstream::out | std::ofstream::binary);
    if (os.fail())
        throw error("htkfile::write(): Open failed");

    // Write the header
    os.write((char*)&nSamples, 4);
    os.write((char*)&sampPeriod, 4);
    os.write((char*)&sampSize, 2);
    os.write((char*)&parmKind, 2);

    // And the rest
    int n = nSamples * sampSize;
    os.write(iVar.ptr<char>(), n);
    if (os.fail())
        throw error("htkfile::write(): Write failed");
}
