/*
 * Copyright 2013 by Idiap Research Institute, http://www.idiap.ch
 *
 * See the file COPYING for the licence associated with this software.
 *
 * Author(s):
 *   Phil Garner, December 2013
 */

#include <lube.h>
#include "ssp.h"
#include "ar.h"

using namespace std;
using namespace ssp;

int main(int argc, char** argv)
{
    var arg(argc, argv);
    cout << "Args: " << arg << endl;

    PCM pcm;
    var a = pcm.read("arctic_a0001.wav");
    std::cout << "a is size: " << a.shape() << std::endl;

    Frame frame(4, 2);
    var f = frame(a);
    cout << "Framed to: " << f.shape() << endl;

    var fv = f.view({4, 4}, 0);
    cout << "fv: " << fv << endl;

    var w = hamming(5);
    w.pop();
    cout << "w: " << w << endl;

    cout << f.atypeStr() << " " << w.atypeStr() << endl;
    f *= w;
    cout << "fv: " << fv << endl;

    lube::DFT dft(4);
    var A = dft(f);
    //var A = lube::abs(F);
    var av = A.view({1,3}, 0);
    cout << "A(0): " << av << endl;

    var B = lube::pow(A,2);
    var bv = B.view({1,3}, 0);
    cout << "B(0): " << bv << endl;

    lube::DFT idft(4, true);
    var C = idft(B);
    var cv = C.view({1,4}, 0);
    cout << "C(0): " << cv << endl;

    Autocorrelation ac(4);
    var AC = ac(f);
    var acv = AC.view({1,4}, 0);
    cout << "ac(0): " << acv << endl;

    A.log();
    cout << "log A(0): " << av << endl;

    var myac;
    myac = 9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.0,1.0,0.0;
    cout << "myac: " << myac << endl;

    Levinson ar(4);
    var myar = ar(myac);
    cout << "myar: " << myar << endl;
    cout << "myac: " << myac << endl;

    Gain gg(4);
    var mygg = gg(myac, myar);
    cout << "mygg: " << mygg << endl;

    var m;
    m = 1.0f, 2, 3, 4,
        1, 2, 3, 4,
        1, 2, 3, 4,
        1, 2, 3, 4,
        1, 2, 3, 4;
    m = m.view({5,4});
    cout << "OLA: " << m << endl;
    OverlapAdd ola;
    var n = ola(m);
    cout << "OLA: " << n << endl;

#if 0
octave:3> a=[1,0.25,0.75]
octave:4> b=[0.5,1.5,2.5]
octave:5> x = [1 2 3 4 5 4 3]
octave:6> filter(b,a,x)
ans =
    0.50000    2.37500    6.03125    8.21094    9.42383   10.98584   10.18567
#endif

    // Make sure these are all asymmetric
    var fa;
    fa = 1.0f, 0.25, 0.75;
    var fb;
    fb = 0.5f, 1.5, 2.5;
    cout << "Numer: " << fb << endl;
    cout << "Denom: " << fa << endl;
    var fx;
    fx = 1.0f, 2, 3, 4, 5, 4, 3;
    Filter ff(fb, fa);
    cout << "Filter: " << ff(fx) << endl;

    var rand = normal(7);
    cout << "Normal: " << rand << endl;

    // LSPs
    // ar: [ 3.63331397 -5.27474314  3.59375591 -0.97826939]
    // ls: [ 0.33440076  0.36424912  0.48876133  0.51274707]
    // ar: [ 3.63331397 -5.27474314  3.59375591 -0.97826939]
    ToLSP tolsp(4);
    FromLSP frlsp(4);
    var pyar;
    pyar = 1.0f, -3.63331397, 5.27474314, -3.59375591, 0.97826939;
    cout << "AR:  " << pyar << endl;
    var lsp = tolsp(pyar);
    cout << "LSP: " << lsp << endl;
    var psl = frlsp(lsp);
    cout << "AR:  " << psl << endl;

    // HTK files
    file htk("htk");
    htk.write("test.htk", fv);
    var fvl = htk.read("test.htk");
    cout << "Read file: " << fvl << endl;
    
    // Done
    return 0;
}
