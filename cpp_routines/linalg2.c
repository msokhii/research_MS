
// Copyright Michael Monagan 2019-2020
// Compile with  gcc -O3 -shared -o linalg2.so -fPIC linalg2.c
// gcc -O3 -c linalg2.c

#define LONG long long int
#define ULONG unsigned long long int
#include <stdio.h>
#include <stdlib.h>
#include "int128g.c"
#include "zmod.c"


/******************************************************************************************/
/*       Zp arithmetic                                                                    */
/******************************************************************************************/
#define add64s linalgadd64s
#define sub64s linalgsub64s
#define neg64s linalgneg64s
#define mul64s linalgmul64s
#define inv64s linalginv64s

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

/********************************************************************************/
/* Linear algebra routines                                                      */
/********************************************************************************/

LONG * matrix64s( int n )
{   LONG *A; LONG N;
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

void matprint64s( LONG *A, int n, int m )
{   int i;
    printf("[");
    for( i=0; i<n; i++ ) { if( i ) printf(" "); vecprint64s(A+i*m,m); printf("\n"); }
    printf("];\n");
    return;
}

LONG det64s( LONG * A, int n, LONG p ) {
// M = [[a,b,c],[d,e,f],[1,2,3],[4,5,6]]  M[i][j]
// A = [a,b,c,d,e,f,1,2,3,4,5,6]          A[m*i+j]   m = #columns 

    int i,j,k;
    LONG d,t,*rowi,*rowk;
    recint P;
    P = recip1(p);
    d = 1;
    for( k=0; k<n; k++ ) {
        //if( k>0 && k%100 == 0 ) printf("elimination at row %d\n",k);
        for( i=k; i<n && A[n*i+k]==0; i++ ); // look for non-zero pivot
        if( i>=n ) { d = 0; break; }
        if( i!=k ) { // interchange row k with row i
             for( j=k; j<n; j++ ) { t = A[k*n+j]; A[k*n+j] = A[i*n+j]; A[i*n+j] = t; }
             d = neg64s(d,p);
        };
        rowk = A+k*n;
        d = mulrec64(d,rowk[k],P); // d = mulrec64(d,A[k*n+k],P);
        //if( A[k*n+k]==0 ) { printf("division by 0\n"); return -1; }
        if( rowk[k]==0 ) { printf("division by 0\n"); return -1; }
        t = inv64s(rowk[k],p); // t = inv64s(A[k*n+k],p);
        rowk[k] = 1; // A[n*k+k] = 1;
        //for( j=k+1; j<n; j++ ) { A[n*k+j] = mulrec64(A[n*k+j],t,P); }
        for( j=k+1; j<n; j++ ) { rowk[j] = mulrec64(rowk[j],t,P); }
        for( i=k+1; i<n; i++ ) { // row i
            rowi = A+i*n;
            if( rowi[k]!=0 ) { // if( A[i*n+k]!=0 )
                for( j=k+1; j<n; j++ ) 
                     // A[i,j] = A[i,j] - A[i,k] A[k,j] mod p
                     // A[i*n+j] = sub64s(A[i*n+j],mulrec64(A[i*n+k],A[k*n+j],P), p );
                     rowi[j] = sub64s(rowi[j],mulrec64(rowi[k],rowk[j],P), p );
                rowi[k] = 0; // A[i*n+k] = 0;
            }
        }
    }
    return d;
}


int minor64s( LONG * A, int n, int m, int *rows, int *cols, LONG p ) {
// Input A[n,m], rows[n], cols[m] 
// Output  r = rank A and rows and cols s.t. rank(B)=r where B=A[rows[0..r-1],cols[0..r-1]]
// Example  A = [[1,1,1],[2,2,2],[1,2,3],[1,0,-1]]
// r = rank(A) = 2 rows = [1,3,-,-], cols = [1,2,-,-] ==> M = [[1,1],[1,2]]
    int i,j,k,r,c,rank;
    LONG t,m64;
    recint P;
    P = recip1(p);
    for( i=0; i<n; i++ ) rows[i] = i;
    rank = 0;
    m64 = m; // so the matrix can have more than 2^31 entries
    r = 0; // current row
    c = 0; // current column
    while ( r<n && c<m ) {
        for( i=r; i<n && A[rows[i]*m64+c]==0; i++ ); // look for non-zero pivot
        if( i>=n ) { c++; continue; } // move to next column
        if( i!=r ) { j = rows[i]; rows[i] = rows[r]; rows[r] = j; }; // row interchange
        cols[rank] = c;
        rank ++;
        k = rows[r]; // pivot row
        t = inv64s(A[k*m64+c],p);
        A[k*m64+c] = 1; // make the pivot row have a pivot = 1
        for( j=c+1; j<m; j++ ) A[k*m64+j] = mulrec64(A[k*m64+j],t,P);
        for( i=r+1; i<n; i++ ) {
            if( A[rows[i]*m64+c]==0 ) continue;
            for( j=c+1; j<m; j++ ) {
                if( A[k*m64+j]==0 ) continue;
                t = mulrec64(A[rows[i]*m64+c],A[k*m64+j],P);
                A[rows[i]*m64+j] = sub64s(A[rows[i]*m64+j],t,p);
            }
            A[rows[i]*m64+c] = 0;
        }
        r++; // next row
        c++; // next column
//printf("r=%d  c=%d\n", r, c);
//printf("A=\n"); matprint64s( A, n, m );
    }
    return rank;
}

LONG dotprod64s( LONG *w, LONG *y, int n, LONG p )
{   // w.y mod p
    int i; LONG s; ULONG z[2];
    //for( i=0,s=0; i<n; i++ ) s = add64s(s,mulrec64(w[i],y[i],P),p);
    z[0] = z[1] = 0;
    for( i=0; i<n-1; ) {
        zfma(z,w[i],y[i]); i++;
        zfma(z,w[i],y[i]); i++;
        if( z[1]>=p ) z[1] -= p;
    }
    if( i==n-1 ) zfma(z,w[i],y[i]);
    zmod(z,p);
    s = z[0];
    return s;
}

void matvecmul64s( LONG *A, int n, LONG *x, LONG *y, LONG p ) {
    int i;
    for( i=0; i<n; i++ ) y[i] = dotprod64s(A+i*n,x,n,p);
    return;
}

LONG matinv64s( LONG * A, int n, LONG * B, LONG p ) {
// Compute A^(-1) in B.  This destroys A.
    int i,j,k;
    LONG d,t;
    recint P;
    P = recip1(p);
    d = 1;
    for( i=0; i<n; i++ ) for( j=0; j<n; j++ ) B[i*n+j] = 0;
    for( i=0; i<n; i++ ) B[i*n+i] = 1;
    for( k=0; k<n; k++ ) {
        for( i=k; i<n && A[n*i+k]==0; i++ ); // look for non-zero pivot
        if( i>=n ) return 0;
        if( i!=k ) { // interchange row k with row i
             for( j=k; j<n; j++ ) { t = A[k*n+j]; A[k*n+j] = A[i*n+j]; A[i*n+j] = t; }
             for( j=0; j<n; j++ ) { t = B[k*n+j]; B[k*n+j] = B[i*n+j]; B[i*n+j] = t; }
             d = neg64s(d,p);
        };
        d = mulrec64(d,A[k*n+k],P);
        t = inv64s(A[k*n+k],p);
        A[n*k+k] = 1;
        for( j=k+1; j<n; j++ ) { A[n*k+j] = mulrec64(A[n*k+j],t,P); }
        for( j=0;   j<n; j++ ) { B[n*k+j] = mulrec64(B[n*k+j],t,P); }
        for( i=0; i<n; i++ ) { // row i
            if( i==k ) continue;
            if( A[i*n+k]!=0 ) {
                for( j=k+1; j<n; j++ ) 
                     if( A[k*n+j] ) A[i*n+j] = sub64s(A[i*n+j],mulrec64(A[i*n+k],A[k*n+j],P), p );
                for( j=0; j<n; j++ )
                     if( B[k*n+j] ) B[i*n+j] = sub64s(B[i*n+j],mulrec64(A[i*n+k],B[k*n+j],P), p );
            }
            A[i*n+k] = 0;
        }
    }
    return d;
}


LONG rref( LONG *B, int n, int m, LONG p ) {
// Put B in reduced row Echelon form and return rank(B)
// The code assumes 0 <= B[i,j] < p
LONG t,det;
int c,r,i,j;
recint P;
   //printf("rref: n=%d m=%d\n",n,m);
   P = recip1(p);
   det = 1;
   for( c=0,r=0; c<m && r<n; c++ ) {
      // Search for a pivot element
      for( i=r; i<n && B[m*i+c]==0; i++ );
      if( i==n ) { det = 0; continue; }
      if( i!=r ) { // interchange row i with row r
          det = neg64s(det,p);
          for( j=c; j<m; j++ ) { t = B[i*m+j]; B[i*m+j] = B[r*m+j]; B[r*m+j] = t; }
      }
      det = mulrec64(det,B[r*m+c],P);
      t = inv64s(B[r*m+c],p);
      for( j=c+1; j<m; j++ ) B[r*m+j] = mulrec64(t,B[r*m+j],P);
      B[r*m+c] = 1;
      for( i=0; i<n; i++ ) {
         if( i==r || B[i*m+c]==0 ) continue;
         for( j=c+1; j<m; j++ )
            B[i*m+j] = sub64s(B[i*m+j],mulrec64(B[i*m+c],B[r*m+j],P),p);
         B[i*m+c] = 0;
      }
      r++;  // go to next row
   }
   printf("det(B)=%lld\n",det);
   return(r); // r = rank(B)
}

/*
int main() {
LONG p, A[12] = {7,3,6,10,2,3,4,5,8,9,0,2};
int i,j,r,n,m;
n = 3;
m = 4;
p = 11;
printf("A = \n"); matprint64s(A,n,m); printf("\n");
r = rref(A,n,m,p);
printf("rank(A)=%d\n",r);
printf("A =\n"); matprint64s(A,n,m); printf("\n");
return 1;
}
*/
