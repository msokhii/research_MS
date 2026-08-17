#include <iostream>
#include <vector>
#include <stdexcept>
#include "polyMath.h"
#include "integerMath.h"

using namespace std;

void vSOLVER64(const vector<LONG>& m,
               const vector<LONG>& y,
               int n,
               vector<LONG>& a,
               vector<LONG>& M,
               int shift,
               LONG p)
{
    if ((int)m.size() < n || (int)y.size() < n) {
        throw invalid_argument("m and y must have at least n entries");
    }

    if (n <= 0) {
        a.clear();
        M.clear();
        return;
    }

    auto normp = [p](LONG x) -> LONG {
        x %= p;
        if (x < 0) x += p;
        return x;
    };

    // Normalize inputs mod p
    vector<LONG> mm(n), yy(n);
    for (int i = 0; i < n; ++i) {
        mm[i] = normp(m[i]);
        yy[i] = normp(y[i]);
    }

    // Distinctness check
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (mm[i] == mm[j]) {
                a.assign(n, 0);
                M.assign(n + 1, 0);
                return;
            }
        }
    }

    a.assign(n, 0);
    M.assign(n + 1, 0);

    // Build Prod(x) = product_i (x - mm[i])
    vector<LONG> Prod(1, 1);   // degree 0
    int degProd = 0;

    vector<LONG> A(2, 0);      // linear factor x - mm[i]
    A[1] = 1;

    for (int i = 0; i < n; ++i) {
        A[0] = neg64s(mm[i], p);
        degProd = pMULIP64VANDER(A, Prod, 1, degProd, p);
    }

    // Return full product polynomial
    for (int i = 0; i <= degProd && i <= n; ++i) {
        M[i] = Prod[i];
    }

    // Sanity: each node is a root of Prod
    for (int j = 0; j < n; ++j) {
        LONG chk = evalHORN64(Prod, mm[j], p);
        if (chk != 0) {
            a.assign(n, 0);
            return;
        }
    }

    // Solve for coefficients
    for (int j = 0; j < n; ++j) {
        // Divide Prod by (x - mm[j])
        vector<LONG> T = Prod;
        vector<LONG> Lin(2, 0);
        Lin[0] = neg64s(mm[j], p);
        Lin[1] = 1;

        int degR = polDIVIP64(T.data(), Lin.data(), degProd, 1, p);
        if (degR != -1) {
            a.assign(n, 0);
            return;
        }

        // Extract quotient Q from T[1..degProd]
        vector<LONG> Q(n, 0);
        for (int i = 0; i < n; ++i) {
            Q[i] = T[i + 1];
        }

        // u = Q(mm[j])
        LONG u = evalHORN64(Q, mm[j], p);
        if (u == 0) {
            a.assign(n, 0);
            return;
        }
        u = modinv64b(u, p);

        // s = sum_i Q[i] * yy[i]
        LONG s = 0;
        for (int i = 0; i < n; ++i) {
            s = add64b(s, mul64bASM(Q[i], yy[i], p), p);
        }

        s = mul64bASM(u, s, p);

        if (shift != 0) {
            LONG invmj = modinv64b(mm[j], p);
            LONG pw = powmod64s(invmj, shift, p);
            s = mul64bASM(pw, s, p);
        }

        a[j] = s;
    }
}