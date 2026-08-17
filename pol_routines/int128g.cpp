#include "int128g.hpp"

// extern long long GLOBALMUL;
volatile int __int128zero = 0;
volatile int __int128junk = 0;
INT64old CNTR = 0;

/* ------------------- functions ------------------- */

void add128(UINT64old *x, UINT64old *y)
{
    ADD22(x[0], x[1], y[0], y[1]);
}

void sub128(UINT64old *x, UINT64old *y)
{
    SUB22(x[0], x[1], y[0], y[1]);
}

void mul128(UINT64old *x, UINT64old a, UINT64old b)
{
    MUL211(x[0], x[1], a, b);
}

void fma128(UINT64old *x, UINT64old a, UINT64old b)
{
    FMA211(x[0], x[1], a, b);
}

void fms128(UINT64old *x, UINT64old a, UINT64old b)
{
    FMS211(x[0], x[1], a, b);
}

void add256(UINT64old *x, UINT64old *y)
{
    ADD44(x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]);
}

void sub256(UINT64old *x, UINT64old *y)
{
    SUB44(x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]);
}

void mul256(UINT64old *x, UINT64old *y, UINT64old *z)
{
    MUL422(x[0], x[1], x[2], x[3], y[0], y[1], z[0], z[1]);
}

void fma256(UINT64old *x, UINT64old *y, UINT64old *z)
{
    FMA422(x[0], x[1], x[2], x[3], y[0], y[1], z[0], z[1]);
}

void fms256(UINT64old *x, UINT64old *y, UINT64old *z)
{
    FMS422(x[0], x[1], x[2], x[3], y[0], y[1], z[0], z[1]);
}

/* ------------------- division ------------------- */

/* 128/64 bit division from Hacker's Delight */
static void _div21(UINT64old u0, UINT64old u1, UINT64old v, UINT64old *q, UINT64old *r)
{
#ifdef DIV21H
    DIV21H(u0, u1, v);
    *q = u0;
    *r = u1;
#else
    UINT64old b, un1, un0, vn1, vn0, q1, q0, un32, un21, un10, rhat;
    INT64old s;

    b = (UINT64old)1 << 32;

    /* quotient overflow case */
    if (u1 > v) DIV_EXCEPTION();
    for (s = 0; !(v >> 63); s++) {
        v = v << 1;
    }
    vn1  = v >> 32;
    vn0  = v & 0xFFFFFFFFULL;
    un32 = (u1 << s) | ((u0 >> (64 - s)) & ((UINT64old)((-s) >> 63)));
    un10 = u0 << s;
    un1  = un10 >> 32;
    un0  = un10 & 0xFFFFFFFFULL;

    q1   = un32 / vn1;
    rhat = un32 - q1 * vn1;
again1:
    if (q1 >= b || q1 * vn0 > b * rhat + un1) {
        q1 = q1 - 1;
        rhat = rhat + vn1;
        if (rhat < b) goto again1;
    }
    un21 = un32 * b + un1 - q1 * v;
    q0   = un21 / vn1;
    rhat = un21 - q0 * vn1;
again2:
    if (q0 >= b || q0 * vn0 > b * rhat + un0) {
        q0 = q0 - 1;
        rhat = rhat + vn1;
        if (rhat < b) goto again2;
    }
    *r = (un21 * b + un0 - q0 * v) >> s;
    *q = q1 * b + q0;
#endif
}

/* [x0:x1] = [x / a : x % a] */
void div128(UINT64old *x, UINT64old a)
{
    _div21(x[0], x[1], a, &x[0], &x[1]);
}

/* return a*b mod p */
UINT64old mulmod64(UINT64old a, UINT64old b, UINT64old p)
{
    UINT64old x0, x1;
    MUL211(x0, x1, a, b);
    _div21(x0, x1, p, &x0, &x1);
    return x1;
}

/* multiword division from Hacker's Delight */
/* requires that n >= m and v[n-1] non-zero */
static void _divmnu(UINT32old *u, UINT32old *v, UINT32old *q, UINT32old *r, INT32old m, INT32old n)
{
    UINT64old b = 4294967296ULL;
    UINT32old un[m + 1], vn[n];
    UINT64old qhat, rhat, p;
    INT64old t, k;
    INT32old s, i, j;

    while (m && u[m - 1] == 0) m--;
    while (n && v[n - 1] == 0) n--;
    if (n == 0) DIV_EXCEPTION();
    if (n == 1) {
        /* single digit divisor */
        for (k = 0, j = m - 1; j >= 0; j--) {
            q[j] = (k * b + u[j]) / v[0];
            k    = (k * b + u[j]) - (UINT64old)q[j] * v[0];
        }
        r[0] = (UINT32old)k;
        return;
    }

    /* normalize divisor to set leading bit */
    for (s = 0, t = v[n - 1]; !(t >> 31); t <<= 1, s++) {}
    for (i = n - 1; i > 0; i--) {
        vn[i] = (v[i] << s) | ((UINT64old)v[i - 1] >> (32 - s));
    }
    vn[0] = v[0] << s;

    /* shift up dividend */
    un[m] = (UINT64old)u[m - 1] >> (32 - s);
    for (i = m - 1; i > 0; i--) {
        un[i] = (u[i] << s) | ((UINT64old)u[i - 1] >> (32 - s));
    }
    un[0] = u[0] << s;

    /* multiply and subtract */
    for (j = m - n; j >= 0; j--) {
        /* compute next quotient */
        p    = un[j + n] * b + un[j + n - 1];
        qhat = p / vn[n - 1];
        rhat = p % vn[n - 1];
again:
        if (qhat >= b || qhat * vn[n - 2] > b * rhat + un[j + n - 2]) {
            qhat--;
            rhat += vn[n - 1];
            if (rhat < b) goto again;
        }
        for (k = 0, i = 0; i < n; i++) {
            p      = qhat * vn[i];
            t      = (INT64old)un[i + j] - k - (INT64old)(p & 0xFFFFFFFFULL);
            un[i + j] = (UINT32old)t;
            k      = (p >> 32) - (UINT64old)(t >> 32);
        }
        t       = (INT64old)un[j + n] - k;
        un[j + n] = (UINT32old)t;
        q[j]    = (UINT32old)qhat;
        if (t < 0) {
            q[j] = q[j] - 1;
            for (k = 0, i = 0; i < n; i++) {
                t       = (UINT64old)un[i + j] + vn[i] + k;
                un[i + j] = (UINT32old)t;
                k       = (UINT64old)t >> 32;
            }
            un[j + n] += (UINT32old)k;
        }
    }

    /* construct remainder */
    for (i = 0; i < n - 1; i++) {
        r[i] = (un[i] >> s) | ((UINT64old)un[i + 1] << (32 - s));
    }
    r[n - 1] = un[n - 1] >> s;
}

/* 256/128 bit long division */
/* quotient/remainder of x/y */
static void _div42(UINT64old *x, UINT64old *y, UINT64old *q, UINT64old *r)
{
    UINT32old R[10] = {0};
    UINT32old Q[10] = {0};

    /* quotient overflow */
    if (x[3] > y[1] || (x[3] == 0 && y[1] == 0 && x[2] > y[0])) {
        DIV_EXCEPTION();
    }
    /* divisor greater than dividend */
    if (x[3] == 0 && x[2] == 0 && (x[1] < y[1] || (x[1] == y[1] && x[0] < y[0]))) {
        q[0] = 0;
        q[1] = 0;
        r[0] = x[0];
        r[1] = x[1];
        return;
    }
    _divmnu((UINT32old *)x, (UINT32old *)y, Q, R, 8, 4);
    r[0] = ((UINT64old)R[1] << 32) | R[0];
    r[1] = ((UINT64old)R[3] << 32) | R[2];
    q[0] = ((UINT64old)Q[1] << 32) | Q[0];
    q[1] = ((UINT64old)Q[3] << 32) | Q[2];
}

/* [x01:x23] = [x / a : x % a] */
void div256(UINT64old *x, UINT64old *y)
{
    UINT64old q[2], r[2];
    _div42(x, y, q, r);
    x[0] = q[0];
    x[1] = q[1];
    x[2] = r[0];
    x[3] = r[1];
}

/* [x0:x1] = [y0:y1] * [z0:z1] mod [p0:p1] */
void mulmod128(UINT64old *x, UINT64old *y, UINT64old *z, UINT64old *p)
{
    UINT64old t[4], q[2], r[2];
    MUL422(t[0], t[1], t[2], t[3], y[0], y[1], z[0], z[1]);
    _div42(t, p, q, r);
    x[0] = r[0];
    x[1] = r[1];
}

/* --------------- inverse division --------------- */

/* for 2/1 division */
static UINT64old _reciprocal21(UINT64old p)
{
    UINT64old u0, u1, q, r;
    u1 = -(p + 1);
    u0 = (UINT64old)-1;
    _div21(u0, u1, p, &q, &r);
    return q;
}

/* for 3/2 division */
static UINT64old _reciprocal32(UINT64old d0, UINT64old d1)
{
    UINT64old t0, t1, p, v;
    v = _reciprocal21(d1);
    p = d1 * v + d0;
    if (p < d0) {
        v--;
        if (p >= d1) {
            v--;
            p -= d1;
        }
        p -= d1;
    }
    MUL211(t0, t1, v, d0);
    p += t1;
    if (p < t1) {
        v--;
        if (p >= d1 || (p == d1 && t0 >= d0)) {
            v--;
        }
    }
    return v;
}

/* 64-bit reciprocal */
recint recip1(UINT64old p)
{
    UINT64old s;
    recint x;
    if (p == 0) DIV_EXCEPTION();
    for (s = 0; !(p >> 63); p <<= 1, s++) {}
    x.s  = s;
    x.v  = _reciprocal21(p);
    x.d0 = p;
    x.d1 = p;
    return x;
}

/* 128-bit reciprocal */
recint recip2(UINT64old *p)
{
    UINT64old p0, p1, t, s;
    recint x;
    p0 = p[0];
    p1 = p[1];
    if (p1 == 0) DIV_EXCEPTION();
    for (s = 0, t = p1; !(t >> 63); t <<= 1, s++) {}
    if (s) {
        p1 = (p1 << s) | (p0 >> (64 - s));
        p0 = (p0 << s);
    }
    x.s  = s;
    x.v  = _reciprocal32(p0, p1);
    x.d0 = p0;
    x.d1 = p1;
    return x;
}

/* 128/64 bit division */
static void _rec21(UINT64old u0, UINT64old u1, UINT64old s, UINT64old v, UINT64old d, UINT64old *q, UINT64old *r)
{
    UINT64old q0, q1, r0;
    if (s) {
        u1 = (u1 << s) | (u0 >> (64 - s));
        u0 = (u0 << s);
    }
    MUL211(q0, q1, v, u1);
    ADD22(q0, q1, u0, u1 + 1);
    r0 = u0 - q1 * d;
    if (r0 > q0) {
        q1--;
        r0 += d;
    }
    if (r0 >= d) {
        q1++;
        r0 -= d;
    }
    q[0] = q1;
    r[0] = r0 >> s;
}

/* 128/64 bit remainder, [u0:u1] already shifted up */
static UINT64old _rec21r(UINT64old u0, UINT64old u1, UINT64old v, UINT64old d)
{
    UINT64old q0, q1, r0;
    MUL211(q0, q1, v, u1);
    ADD22(q0, q1, u0, u1 + 1);
    r0 = u0 - q1 * d;
    if (r0 > q0) {
        r0 += d;
    }
    if (r0 >= d) {
        r0 -= d;
    }
    return r0;
}

/* 192/128 bit division */
static void _rec32(UINT64old u0, UINT64old u1, UINT64old u2, UINT64old s, UINT64old v,
                   UINT64old d0, UINT64old d1, UINT64old *q, UINT64old *r)
{
    UINT64old q0, q1, r0, r1, t0, t1;
    if (s) {
        u2 = (u2 << s) | (u1 >> (64 - s));
        u1 = (u1 << s) | (u0 >> (64 - s));
        u0 = (u0 << s);
    }
    MUL211(q0, q1, v, u2);
    ADD22(q0, q1, u1, u2);
    r0 = u0;
    r1 = u1 - q1 * d1;
    MUL211(t0, t1, d0, q1);
    SUB22(r0, r1, d0, d1);
    SUB22(r0, r1, t0, t1);
    q1++;
    if (r1 >= q0) {
        q1--;
        ADD22(r0, r1, d0, d1);
    }
    if (r1 > d1 || ((r1 == d1) && (r0 >= d0))) {
        q1++;
        SUB22(r0, r1, d0, d1);
    }
    if (s) {
        r0 = (r0 >> s) | (r1 << (64 - s));
        r1 = (r1 >> s);
    }
    q[0] = q1;
    r[0] = r0;
    r[1] = r1;
}

/* 256/128 bit division */
static void _rec42(UINT64old u0, UINT64old u1, UINT64old u2, UINT64old u3, UINT64old s, UINT64old v,
                   UINT64old d0, UINT64old d1, UINT64old *q, UINT64old *r)
{
    UINT64old q0, q1, R[3];
    _rec32(u1, u2, u3, s, v, d0, d1, &q1, R + 1);
    _rec32(u0, R[1], R[2], s, v, d0, d1, &q0, R);
    q[0] = q0;
    q[1] = q1;
    r[0] = R[0];
    r[1] = R[1];
}

/* 128/64 bit division */
void rec128(UINT64old *x, recint v)
{
    _rec21(x[0], x[1], v.s, v.v, v.d0, &x[0], &x[1]);
}

/* 64 x 64 bit mulmod */
UINT64old mulrec64(UINT64old a, UINT64old b, recint v)
{   
    // ++GLOBALMUL;
    UINT64old u0, u1, r;
    b = b << v.s;
    MUL211(u0, u1, a, b);
    r = _rec21r(u0, u1, v.v, v.d0);
    return r >> v.s;
}

/* 256/128 bit division */
void rec256(UINT64old *x, recint v)
{
    UINT64old q[2], r[2];
    _rec42(x[0], x[1], x[2], x[3], v.s, v.v, v.d0, v.d1, q, r);
    x[0] = q[0];
    x[1] = q[1];
    x[2] = r[0];
    x[3] = r[1];
}

/* 128 x 128 bit mulmod */
void mulrec128(UINT64old *x, UINT64old *y, UINT64old *z, recint v)
{
    UINT64old t[4], q[2], r[2];
    MUL422(t[0], t[1], t[2], t[3], y[0], y[1], z[0], z[1]);
    _rec42(t[0], t[1], t[2], t[3], v.s, v.v, v.d0, v.d1, q, r);
    x[0] = r[0];
    x[1] = r[1];
}

/* --------------- 128-bit integers --------------- */
#if INT128_EXTENDED

void mul256i(UINT64old *x, UINT128 a, UINT128 b)
{
    UINT64old A[2], B[2];
    A[0] = *((UINT64old *)(&a) + 0);
    A[1] = *((UINT64old *)(&a) + 1);
    B[0] = *((UINT64old *)(&b) + 0);
    B[1] = *((UINT64old *)(&b) + 1);
    MUL422(x[0], x[1], x[2], x[3], A[0], A[1], B[0], B[1]);
}

void div256i(UINT64old *x, UINT128 a, UINT128 *q, UINT128 *r)
{
    UINT128 q1, r1;
    UINT64old A[2], Q[2], R[2];
    A[0] = *((UINT64old *)(&a) + 0);
    A[1] = *((UINT64old *)(&a) + 1);
    _div42(x, A, Q, R);
    q1 = Q[1];
    q1 = (q1 << 64) | Q[0];
    r1 = R[1];
    r1 = (r1 << 64) | R[0];
    *q = q1;
    *r = r1;
}

UINT128 mulmod128i(UINT128 a, UINT128 b, UINT128 p)
{
    UINT128 r;
    UINT64old A[2], B[2], P[2], Q[2] = {0, 0}, R[2] = {0, 0}, x[4] = {0, 0, 0, 0};
    A[0] = *((UINT64old *)(&a) + 0);
    A[1] = *((UINT64old *)(&a) + 1);
    B[0] = *((UINT64old *)(&b) + 0);
    B[1] = *((UINT64old *)(&b) + 1);
    P[0] = *((UINT64old *)(&p) + 0);
    P[1] = *((UINT64old *)(&p) + 1);
    MUL422(x[0], x[1], x[2], x[3], A[0], A[1], B[0], B[1]);
    _div42(x, P, Q, R);
    r = R[1];
    r = (r << 64) | R[0];
    return r;
}

UINT128 mulrec128i(UINT128 a, UINT128 b, recint v)
{
    UINT128 r;
    UINT64old A[2], B[2], Q[2] = {0, 0}, R[2] = {0, 0}, x[4] = {0, 0, 0, 0};
    A[0] = *((UINT64old *)(&a) + 0);
    A[1] = *((UINT64old *)(&a) + 1);
    B[0] = *((UINT64old *)(&b) + 0);
    B[1] = *((UINT64old *)(&b) + 1);
    MUL422(x[0], x[1], x[2], x[3], A[0], A[1], B[0], B[1]);
    _rec42(x[0], x[1], x[2], x[3], v.s, v.v, v.d0, v.d1, Q, R);
    r = R[1];
    r = (r << 64) | R[0];
    return r;
}

#endif
