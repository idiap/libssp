/*
 * Copyright 2016 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, August 2016
 */

#ifndef WINDOW_H
#define WINDOW_H

#include <ssp/ssp.h>

namespace ssp
{
    class Window : public UnaryFunctor
    {
    public:
        /**
         * Set window type
         *
         * See http://en.wikipedia.org/wiki/Window_function
         */
        virtual void
        set(int iSize, bool iPeriodic=false, float iParam=0.0f) = 0;
        Window(int iSize) : UnaryFunctor(iSize) {};
        explicit operator var () const { return mWindow; };
    protected:
        var mWindow;
    };

    class Hann : public Window
    {
    public:
        Hann(int iSize, bool iPeriodic=false)
            : Window(iSize) { set(iSize, iPeriodic); };
        virtual void set(int iSize, bool iPeriodic=false, float iParam=0.0f);
    };

    class Hamming : public Window
    {
    public:
        Hamming(int iSize, bool iPeriodic=false)
            : Window(iSize) { set(iSize, iPeriodic); };
        virtual void set(int iSize, bool iPeriodic=false, float iParam=0.0f);
    };

    class Nuttall : public Window
    {
    public:
        Nuttall(int iSize, bool iPeriodic=false)
            : Window(iSize) { set(iSize, iPeriodic); };
        virtual void set(int iSize, bool iPeriodic=false, float iParam=0.0f);
    };

    class BlackmanHarris : public Window
    {
    public:
        BlackmanHarris(int iSize, bool iPeriodic=false)
            : Window(iSize) { set(iSize, iPeriodic); };
        virtual void set(int iSize, bool iPeriodic=false, float iParam=0.0f);
    };

    class BlackmanNuttall : public Window
    {
    public:
        BlackmanNuttall(int iSize, bool iPeriodic=false)
            : Window(iSize) { set(iSize, iPeriodic); };
        virtual void set(int iSize, bool iPeriodic=false, float iParam=0.0f);
    };

    class Gaussian : public Window
    {
    public:
        Gaussian(int iSize, bool iPeriodic=false, float iParam=1.0f)
            : Window(iSize) { set(iSize, iPeriodic, iParam); };
        virtual void set(int iSize, bool iPeriodic=false);
        void set(int iSize, bool iPeriodic=false, float iParam=0.0f);
    };
}

#endif // WINDOW_H
