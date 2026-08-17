#pragma once

#define INT128_ASSEMBLY 1
#define INT128_EXTENDED 1

#define INT32old int
#define INT64old long long
#define UINT32old unsigned int
#define UINT64old unsigned long long

typedef struct {
    UINT64old s;
    UINT64old v;
    UINT64old d0;
    UINT64old d1;
} recint;

#if INT128_EXTENDED
    #define INT128  __int128_t
    #define UINT128 __uint128_t
#endif

/* declarations only */
void add128(UINT64old *x, UINT64old *y);
void add256(UINT64old *x, UINT64old *y);
void sub128(UINT64old *x, UINT64old *y);
void sub256(UINT64old *x, UINT64old *y);

void mul128(UINT64old *x, UINT64old a, UINT64old b);
void mul256(UINT64old *x, UINT64old *y, UINT64old *z);

void fma128(UINT64old *x, UINT64old a, UINT64old b);
void fms128(UINT64old *x, UINT64old a, UINT64old b);
void fma256(UINT64old *x, UINT64old *y, UINT64old *z);
void fms256(UINT64old *x, UINT64old *y, UINT64old *z);

void div128(UINT64old *x, UINT64old a);
void div256(UINT64old *x, UINT64old *y);

UINT64old mulmod64(UINT64old a, UINT64old b, UINT64old p);
void mulmod128(UINT64old *x, UINT64old *y, UINT64old *z, UINT64old *p);

recint recip1(UINT64old p);
recint recip2(UINT64old *p);

void rec128(UINT64old *x, recint v);
void rec256(UINT64old *x, recint v);

UINT64old mulrec64(UINT64old a, UINT64old b, recint v);
void mulrec128(UINT64old *x, UINT64old *y, UINT64old *z, recint v);

#if INT128_EXTENDED
void mul256i(UINT64old *x, UINT128 a, UINT128 b);
void div256i(UINT64old *x, UINT128 a, UINT128 *q, UINT128 *r);
UINT128 mulmod128i(UINT128 a, UINT128 b, UINT128 p);
UINT128 mulrec128i(UINT128 a, UINT128 b, recint v);
#endif

/* globals declared, not defined */
extern volatile int __int128zero;
extern volatile int __int128junk;
extern INT64old CNTR;

/* compile-time check */
static_assert(sizeof(UINT64old) == 8, "UINT64old must be 64-bit");

/* exception macro */
#define DIV_EXCEPTION() __int128junk = 128 / __int128zero

/* -------------------- macros -------------------- */

#define INT128_MACROS 0

#if INT128_ASSEMBLY && defined(__GNUC__) && defined(__x86_64__)

    #define __ADD(x,y) " addq " #x ", " #y "\n\t"
    #define __ADC(x,y) " adcq " #x ", " #y "\n\t"
    #define __SUB(x,y) " subq " #x ", " #y "\n\t"
    #define __SBB(x,y) " sbbq " #x ", " #y "\n\t"
    #define __MUL(x)   " mulq " #x "\n\t"
    #define __DIV(x)   " divq " #x "\n\t"

    #define ADD22(x0,x1,y0,y1) __asm__(\
        __ADD(%2, %0) \
        __ADC(%3, %1) \
    : "=r"(x0), "=r"(x1) : "g"(y0), "g"(y1), "0"(x0), "1"(x1) : "cc")

    #define ADD33(x0,x1,x2,y0,y1,y2) __asm__(\
        __ADD(%3, %0) \
        __ADC(%4, %1) \
        __ADC(%5, %2) \
    : "=r"(x0), "=r"(x1), "=r"(x2) : "g"(y0), "g"(y1), "g"(y2), "0"(x0), "1"(x1), "2"(x2) : "cc")

    #define ADD44(x0,x1,x2,x3,y0,y1,y2,y3) __asm__(\
        __ADD(%4, %0) \
        __ADC(%5, %1) \
        __ADC(%6, %2) \
        __ADC(%7, %3) \
    : "=r"(x0), "=r"(x1), "=r"(x2), "=r"(x3) : "g"(y0), "g"(y1), "g"(y2), "g"(y3), "0"(x0), "1"(x1), "2"(x2), "3"(x3) : "cc")

    #define SUB22(x0,x1,y0,y1) __asm__(\
        __SUB(%2, %0) \
        __SBB(%3, %1) \
    : "=r"(x0), "=r"(x1) : "g"(y0), "g"(y1), "0"(x0), "1"(x1) : "cc")

    #define SUB33(x0,x1,x2,y0,y1,y2) __asm__(\
        __SUB(%3, %0) \
        __SBB(%4, %1) \
        __SBB(%5, %2) \
    : "=r"(x0), "=r"(x1), "=r"(x2) : "g"(y0), "g"(y1), "g"(y2), "0"(x0), "1"(x1), "2"(x2) : "cc")

    #define SUB44(x0,x1,x2,x3,y0,y1,y2,y3) __asm__(\
        __SUB(%4, %0) \
        __SBB(%5, %1) \
        __SBB(%6, %2) \
        __SBB(%7, %3) \
    : "=r"(x0), "=r"(x1), "=r"(x2), "=r"(x3) : "g"(y0), "g"(y1), "g"(y2), "g"(y3), "0"(x0), "1"(x1), "2"(x2), "3"(x3) : "cc")

    #define MUL211(x0,x1,a,b) __asm__(\
        __MUL(%3) \
    : "=a"(x0), "=d"(x1) : "0"(a), "r"(b) : "cc")

    #define DIV21H(x0,x1,a) __asm__(\
        __DIV(%2) \
    : "=a"(x0), "=d"(x1) : "r"(a), "0"(x0), "1"(x1) : "cc")

#else

    #define ADD22(x0,x1,y0,y1) do { \
        x0 += y0; \
        x1 += y1 + (x0 < y0); \
    } while (0)

    #define ADD33(x0,x1,x2,y0,y1,y2) do { \
        UINT64old s1; \
        s1 = x1; \
        x0 += y0; \
        x1 += (x0 < y0); \
        x2 += (x1 < s1); \
        x1 += y1; \
        x2 += y2 + (x1 < y1); \
    } while (0)

    #define ADD44(x0,x1,x2,x3,y0,y1,y2,y3) do { \
        UINT64old s1,s2; \
        s1 = x1; s2 = x2; \
        x0 += y0; \
        x1 += (x0 < y0); \
        x2 += (x1 < s1); \
        x3 += (x2 < s2); \
        s2 = x2; \
        x1 += y1; \
        x2 += (x1 < y1); \
        x3 += (x2 < s2); \
        x2 += y2; \
        x3 += y3 + (x2 < y2); \
    } while (0)

    #define SUB22(x0,x1,y0,y1) do { \
        UINT64old s0; \
        s0 = x0; \
        x0 -= y0; \
        x1 -= y1 + (s0 < y0); \
    } while (0)

    #define SUB33(x0,x1,x2,y0,y1,y2) do { \
        UINT64old s0,s1; \
        s0 = x0; s1 = x1; \
        x0 -= y0; \
        x1 -= (s0 < y0); \
        x2 -= (s1 < x1); \
        s1 = x1; \
        x1 -= y1; \
        x2 -= y2 + (s1 < y1); \
    } while (0)

    #define SUB44(x0,x1,x2,x3,y0,y1,y2,y3) do { \
        UINT64old s0,s1,s2; \
        s0 = x0; s1 = x1; s2 = x2; \
        x0 -= y0; \
        x1 -= (s0 < y0); \
        x2 -= (s1 < x1); \
        x3 -= (s2 < x2); \
        s1 = x1; s2 = x2; \
        x1 -= y1; \
        x2 -= (s1 < y1); \
        x3 -= (s2 < x2); \
        s2 = x2; \
        x2 -= y2; \
        x3 -= y3 + (s2 < y2); \
    } while (0)

    #define MUL211(x0,x1,a,b) do { \
        UINT64old a0,a1,b0,b1,t0,t1; \
        a0 = a & 0xFFFFFFFFULL; \
        b0 = b & 0xFFFFFFFFULL; \
        x0 = a0*b0; \
        a1 = a >> 32; \
        b1 = b >> 32; \
        x1 = a1*b1; \
        t0 = a0*b1; \
        t1 = t0 >> 32; \
        t0 = t0 << 32; \
        ADD22(x0,x1,t0,t1); \
        t0 = a1*b0; \
        t1 = t0 >> 32; \
        t0 = t0 << 32; \
        ADD22(x0,x1,t0,t1); \
    } while (0)

#endif

#define FMA211(x0,x1,a,b) do { \
    UINT64old u0,u1; \
    MUL211(u0,u1,a,b); \
    ADD22(x0,x1,u0,u1); \
} while (0)

#define FMS211(x0,x1,a,b) do { \
    UINT64old u0,u1; \
    MUL211(u0,u1,a,b); \
    SUB22(x0,x1,u0,u1); \
} while (0)

#define MUL321(x0,x1,x2,y0,y1,a) do { \
    x2 = 0; \
    MUL211(x0,x1,y0,a); \
    FMA211(x1,x2,y1,a); \
} while (0)

#define FMA321(x0,x1,x2,y0,y1,a) do { \
    UINT64old u2,u3,u4; \
    MUL321(u2,u3,u4,y0,y1,a); \
    ADD33(x0,x1,x2,u2,u3,u4); \
} while (0)

#define FMS321(x0,x1,x2,y0,y1,a) do { \
    UINT64old u2,u3,u4; \
    MUL321(u2,u3,u4,y0,y1,a); \
    SUB33(x0,x1,x2,u2,u3,u4); \
} while (0)

#define MUL422(x0,x1,x2,x3,y0,y1,z0,z1) do { \
    UINT64old t2,t3; \
    MUL211(x0,x1,y0,z0); \
    MUL211(x2,x3,y1,z1); \
    MUL211(t2,t3,y1,z0); \
    ADD33(x1,x2,x3,t2,t3,0); \
    MUL211(t2,t3,y0,z1); \
    ADD33(x1,x2,x3,t2,t3,0); \
} while (0)

#define FMA422(x0,x1,x2,x3,y0,y1,z0,z1) do { \
    UINT64old u5,u6,u7,u8; \
    MUL422(u5,u6,u7,u8,y0,y1,z0,z1); \
    ADD44(x0,x1,x2,x3,u5,u6,u7,u8); \
} while (0)

#define FMS422(x0,x1,x2,x3,y0,y1,z0,z1) do { \
    UINT64old u5,u6,u7,u8; \
    MUL422(u5,u6,u7,u8,y0,y1,z0,z1); \
    SUB44(x0,x1,x2,x3,u5,u6,u7,u8); \
} while (0)

#define MUL431(x0,x1,x2,x3,y0,y1,y2,a) do { \
    UINT64old t4,t5; \
    MUL211(x0,x1,y0,a); \
    MUL211(t4,t5,y1,a); \
    MUL211(x2,x3,y2,a); \
    ADD33(x1,x2,x3,t4,t5,0); \
} while (0)

#define FMA431(x0,x1,x2,x3,y0,y1,y2,a) do { \
    UINT64old u5,u6,u7,u8; \
    MUL431(u5,u6,u7,u8,y0,y1,y2,a); \
    ADD44(x0,x1,x2,x3,u5,u6,u7,u8); \
} while (0)

#define FMS431(x0,x1,x2,x3,y0,y1,y2,a) do { \
    UINT64old u5,u6,u7,u8; \
    MUL431(u5,u6,u7,u8,y0,y1,y2,a); \
    SUB44(x0,x1,x2,x3,u5,u6,u7,u8); \
} while (0)