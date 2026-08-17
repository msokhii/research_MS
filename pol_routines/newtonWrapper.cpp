#include <vector>
#include <cstdint>
#include <algorithm>
#include "integerMath.h"
#include "int128g.hpp"
#include "interpAlgo.h"

using namespace std;
using LONG = int64_t;

extern "C" int cppInterp(int xLen,
                                    const LONG *xIn,
                                    int yLen,
                                    const LONG *yIn,
                                    const LONG p,
                                    int outLen,
                                    LONG *yOut,
                                    int *degOut)
{
    // basic checks
    if (!xIn || !yIn || !yOut || !degOut) {
        return -1;
    }
    if (xLen <= 0 || yLen <= 0 || outLen <= 0) {
        return -2;
    }
    if (xLen != yLen) {
        return -3;
    }
    if (outLen < yLen) {
        return -4;
    }

    const int n = xLen;

    // initialize outputs
    *degOut = -1;
    for (int i = 0; i < outLen; ++i) {
        yOut[i] = 0;
    }

    // local working copies since kernel overwrites y
    vector<LONG> x(xIn, xIn + n);
    vector<LONG> y(yIn, yIn + n);

    recint P = recip1(p);

    int d = newtonInterpMulRec(x.data(), y.data(), n, p, P);

    if (d < 0) {
        *degOut = -1;
        return d;
    }

    if (d >= outLen) {
        *degOut = -1;
        return -5;
    }

    std::copy_n(y.data(), d + 1, yOut);
    *degOut = d;

    return 0;
}