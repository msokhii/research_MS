#include <vector>
#include"polyMath.h"
#include"integerMath.h"
#include"fastVSolve.h"

extern "C" int vSOLVER64C(const LONG *RR,
                          const LONG *LL,
                          int nn,
                          LONG *aa,
                          LONG *XX,
                          int shift,
                          LONG pp)
{
    if (!RR || !LL || !aa || !XX) return -1;
    if (nn < 0) return -2;

    vector<LONG> m(RR, RR + nn);
    vector<LONG> y(LL, LL + nn);
    vector<LONG> a;
    vector<LONG> M;

    try {
        vSOLVER64(m, y, nn, a, M, shift, pp);
    } catch (...) {
        return -3;
    }

    for (int i = 0; i < nn; ++i) {
        aa[i] = (i < (int)a.size() ? a[i] : 0);
    }
    for (int i = 0; i <= nn; ++i) {
        XX[i] = (i < (int)M.size() ? M[i] : 0);
    }

    return 0;
}