/*
 * Copyright 2015 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, February 2015
 */

#ifndef PITCH_H
#define PITCH_H

#include "ssp.h"
#include "ar.h"

namespace ssp
{
    /**
     * Pitch is only a functor in order to take advantage of the allocator
     * (which PitchHNR would do anyway); otherwise it just implements scalar()
     * and keeps track of some variables.  The contents of scalar() could just
     * be written as a function.
     */
    class Pitch : public ssp::UnaryFunctor
    {
    public:
        Pitch(PCM* iPCM, var iLo = 40.0f, var iHi = 500.0f);
    protected:
        void scalar(const var& iVar, var& oVar) const;
        PCM* mPCM;
        var mLo;
        var mHi;
    };

    var excitation(var iPitch, var iHNR, PCM* iPCM);
};

#endif // PITCH_H
