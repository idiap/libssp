/*
 * Copyright 2016 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, January 2016
 */

#ifndef WARP_H
#define WARP_H

namespace ssp
{
    float hzToERB(float iHz)
    {
        float erb = 24.7f * (4.37e-3*iHz + 1.0f);
        return erb;
    }

    float erbToHz(float iERB)
    {
        float hz = (iERB/24.7f-1.0f) / 4.37e-3;
        return hz;
    }

    float hzToERBRate(float iHz)
    {
        float rate = (1000.0f/107.94f)*std::log(437.0e-3f*iHz+1);
        return rate;
    }

    float erbRateToHz(float iRate)
    {
        float hz = (std::exp(107.94e-3*iRate) - 1.0f) / 437e-3f;
        return hz;
    }

    /** Convert a value in Hz to the mel scale. */
    float hzToMel(float iHz)
    {
        return 2595.0f * std::log10(1.0f + iHz / 700.0f);
    }

    /** Convert a value from the mel scale to Hz. */
    float melToHz(float iMel)
    {
        return 700.0 * (std::pow(10, iMel / 2595.0f) - 1.0f);
    }
}

#endif // WARP_H
