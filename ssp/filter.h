/*
 * Copyright 2016 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, April 2016
 */

#ifndef FILTER_H
#define FILTER_H

namespace ssp
{
    namespace core
    {
        /**
         * Filter class
         */
        class Filter
        {
        public:
            Filter();
            Filter(int iNNumer, float* iNumer, int iNDenom, float* iDenom);
            ~Filter();
            void set(int iNNumer, float* iNumer, int iNDenom, float* iDenom);
            void operator()(int iNSamples, float* iSample, float* oSample)
                const;
            float operator()(float iSample, float* iState) const;
        private:
            int mNNumer;
            float* mNumer;
            int mNDenom;
            float* mDenom;
            int mNStates;
        };
    }
}

#endif // FILTER_H
