
//  gcc -O3 -shared -o VSolve1.so -fPIC VSolve1.c polyalg2.c

#include <stdio.h>
#include <time.h>
#define LONG long long int
#define ULONG unsigned long long int


#define UINT64 unsigned long long

typedef struct {
        UINT64 s;       /* shift */
        UINT64 v;       /* reciprocal of d */
        UINT64 d0;      /* divisor shifted up */
        UINT64 d1;
} recint;

recint recip1(UINT64  p);
UINT64 mulrec64(UINT64  a, UINT64  b, recint  v);


/******************************************************************************************/
/*  Zp utilities                                                                          */
/******************************************************************************************/



#define add64s VSadd64s
#define sub64s VSsub64s
#define neg64s VSneg64s
#define mul64s VSmul64s
LONG add64s(LONG a, LONG b, LONG p) { LONG t; t = (a-p)+b; t += (t>>63) & p; return t; }
LONG sub64s(LONG a, LONG b, LONG p) { LONG t; t = a-b; t += (t>>63) & p; return t; }
LONG neg64s(LONG a, LONG p) { return (a==0) ? 0 : p-a; }
LONG mul64s(LONG a, LONG b, LONG p) {
        LONG q, r;
        __asm__ __volatile__(           \
        "       mulq    %%rdx           \n\t" \
        "       divq    %4              \n\t" \
        : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
        return r;
}

LONG modinv64s(LONG u, LONG p);
LONG powmod64s(LONG u, LONG n, LONG p);

int polmul64s( LONG * A, LONG * B, LONG * C, int da, int db, LONG p);
int poldiv64s( LONG * A, LONG * B, int da, int db, LONG p );
LONG poleval64s(LONG *a, int d, LONG x, LONG p);



void VandermondeSolve64s( LONG *m, LONG *y, int n, LONG *a, LONG *M, int shift, LONG p )
{   // m,y,a are dimension n, M is dimension n+1
    // 
    // m = [m1,m2,...,mn] and v = [v1,v2,...,vn]
    // Solve V a = v for a where
    //
    //     [ 1        1        ... 1        ] [ a1 ]   [ y1 ]
    //     [ m1       m2       ... mn       ] [ a2 ]   [ y2 ]
    // V = [ m1^2     m2^2     ... mn^2     ] [ a3 ] = [ y3 ]
    //     [ ...      ...      ... ...      ] [    ]   [    ]
    //     [ m1^(n-1) m2^(n-1) ... mn^(n-1) ] [ an ]   [ yn ]
    //
    // Let M = (x-m[1])(x-m[2])*...*(x-m[n])
    // Let Qj(x) = L(x)/(x-m[j]) = q1 + q2 x + ... + qn x^(n-1).
    // Let u = Qj(m[j]) and v = [q1,q2,...,qn=1] be the array of coefficients of Qj(x).
    // Then the solution a[j] = v dot y / u for 1 <= j <= n.
    //

    int i,j;
    LONG u,s,A[2];
    clock_t st,et;
    recint P;
    
    P = recip1(p);
    M[0] = 1;
    A[1] = 1;
    // compute M = (x-m[0])(x-m[1])...(x-m[n-1])
st = clock();
    for( i=0; i<n; i++ ) { 
         A[0] = neg64s(m[i],p);
         polmul64s( A, M, M, 1, i, p );
    }
et = clock();
printf("Mul time=%lld ms\n", (et-st)/1000 );
    for( j=0; j<n; j++ ) {
         A[0] = neg64s(m[j],p);
         i = poldiv64s(M,A,n,1,p); // M = M / (x-m[j])
         for( i=0; i<n; i++ ) M[i] = M[i+1]; // move quotient to front of M
         u = poleval64s(M,n-1,m[j],p);
         if( u==0 ) { printf("roots are not distinct\n"); return; }
         u = modinv64s(u,p);
         // compute a[j] = M dot y
         //for( s=0,i=0; i<n; i++ ) s = add64s(s,mul64s(M[i],y[i],p),p);
         for( s=0,i=0; i<n; i++ ) s = add64s(s,mulrec64(M[i],y[i],P),p);
         s = mul64s(u,s,p);
         //s = mulrec64(u,s,P);
         if( shift!=0 ) {
             u = modinv64s(m[j],p);
             u = powmod64s(u,shift,p);
             //s = mulrec64(u,s,P);
             s = mul64s(u,s,p);
         }
         a[j] = s;
         polmul64s(A,M,M,1,n-1,p); // restore M = M x (x-m[j])
    }
et = clock();
printf("Solve time=%lld ms\n", (et-st)/1000 );
    return;
}



