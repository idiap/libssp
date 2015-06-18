/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, May 2015
 */

#ifndef ARCODEC_H
#define ARCODEC_H

#include "ssp.h"

namespace ssp
{
    /**
     * AR codec
     */
    class ARCodec : public Codec
    {
    public:
        ARCodec(PCM* iPCM, bool iOracle=false);
        virtual var encode(var iSignal);
        virtual var decode(var iParams);
        virtual var read(var iFile);
        virtual void write(var iFile, var iParams);
    private:
        bool mOracle;
    };
}

#endif // ARCODEC_H
