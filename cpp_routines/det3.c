
// Copyright Michael Monagan 2019
// Compile with  gcc -O3 -shared -o det3.so -fPIC det3.c

#define LONG long long int
#include <stdio.h>
#include <stdlib.h>
#include "int128g.c"



/******************************************************************************************/

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

/* c^(-1) mod p assuming 0 < c < p < 2^63 */
LONG inv64s( LONG c, LONG p )
{   LONG d,r,q,r1,c1,d1;
    d = p;
    c1 = 1;
    d1 = 0;
    while( d != 0 ) {
        q = c / d;
        r = c - q*d; r1 = c1 - q*d1;
        c = d; c1 = d1;
        d = r; d1 = r1;
    }
    if( c!=1 ) return( 0 );
    if( c1 < 0 ) c1 += p;
    return( c1 );
}

/******************************************************************************************/
/* Linear algebra routines                                                                */
/******************************************************************************************/

LONG * matrix( int n )
{
    LONG *A; LONG N;
    N = n;
    N = n*N;
    N = sizeof(LONG) * N;
    A = (LONG *) malloc(N);
    return A;
}

/* print an array in form [a0,a1,...,an-1] */
void vecprint64s( LONG *A, int n )
{   int i;
    printf("[");
    for( i=0; i<n; i++ ) { printf("%lld",A[i]); if( i<n-1 ) printf(", "); }
    printf("]");
    return;
}

void matprint64s( LONG *A, int n )
{   int i;
    printf("[");
    for( i=0; i<n-1; i++ ) { vecprint64s(A+i*n,n); printf(",\n"); }
    vecprint64s(A+(n-1)*n,n); 
    printf("];\n");
    return;
}


/******************************************************************************************/
LONG det64s( LONG * A, int n, LONG p ) {
    int i,j,k,m;
    LONG d,t,*S;
    recint P;
    P = recip1(p);
    d = 1;
    S = (LONG *) malloc(n*sizeof(LONG));
    for( k=0; k<n; k++ ) {
        //if( k>0 && k%100 == 0 ) printf("elimination at row %d\n",k);
        for( i=k; i<n && A[n*i+k]==0; i++ ); // look for non-zero pivot
        if( i>=n ) { d = 0; break; }
        if( i!=k ) { // interchange row k with row i
             for( j=k; j<n; j++ ) { t = A[k*n+j]; A[k*n+j] = A[i*n+j]; A[i*n+j] = t; }
             d = neg64s(d,p);
        };
        d = mulrec64(d,A[k*n+k],P);
        t = inv64s(A[k*n+k],p);
        if( t==0 ) { printf("division by 0\n"); return -1; }
        A[n*k+k] = 1;
        for( m=0,j=k+1; j<n; j++ )
            if( A[n*k+j] ) { A[n*k+j] = mulrec64(A[n*k+j],t,P); S[m]=j; m++; }
        //printf("n-k=%d  m=%d\n",n-k,m);
        for( i=k+1; i<n; i++ ) { // row i
            if( A[i*n+k]!=0 )
                for( j=0; j<m; j++ ) 
                       A[i*n+S[j]] = sub64s(A[i*n+S[j]],mulrec64(A[i*n+k],A[k*n+S[j]],P), p );
                //for( j=k+1; j<n; j++ ) 
                //     A[i*n+j] = sub64s(A[i*n+j],mulrec64(A[i*n+k],A[k*n+j],P), p );
            A[i*n+k] = 0;
        }
    }
    free(S);
    return d;
}
