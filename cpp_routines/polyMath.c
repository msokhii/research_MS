
#define LONG long long int
#define ULONG unsigned long long int


ULONG seed = 1;
ULONG mult = 6364136223846793003ll;

LONG rand64s(LONG p) {
// output a random number on [0,p) using x[k+1] = mult x[k] mod 2^64 with x[0] = 1 = seed
    LONG x,y;
    extern ULONG seed, mult;
    seed = mult*seed;
    x = seed >> 32;
    seed = mult*seed;
    y = seed >> 32;
    x = (x<<31) | y;
    x = x % p;
    return(x);
}


/******************************************************************************************/
/*  Zp utilities   Assumes inputs a and b are in [0,p)                                    */
/******************************************************************************************/


LONG add64s(LONG a, LONG b, LONG p) { LONG t; t = (a-p)+b; t += (t>>63) & p; return t; }
LONG sub64s(LONG a, LONG b, LONG p) { LONG t; t = a-b; t += (t>>63) & p; return t; }
LONG neg64s(LONG a, LONG p) { return (a==0) ? 0 : p-a; }
LONG mod64s(LONG a, LONG p) { 
    if( a<=-p || a>=p ) a = a%p; 
    if( a<0 ) a += p; 
    return a;
}

LONG mul64s(LONG a, LONG b, LONG p) {
        LONG q, r;
        __asm__ __volatile__(           \
        "       mulq    %%rdx           \n\t" \
        "       divq    %4              \n\t" \
        : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
        return r;
}


/* c^(-1) mod p assuming 0 < c < p < 2^63 */
LONG modinv64s( LONG c, LONG p )
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


/* a^n mod p assuming 0 <= a < p < 2^63 */
LONG powmod64s( LONG a, LONG n, LONG p )
{   LONG r,s;
    a += (a>>63) & p; // protect from bad input
    if( n==0 ) return 1;
    if( n==1 ) return a;
    for( r=1, s=a; n>0; n /= 2 ) { if( n & 1 ) r = mul64s(r,s,p); s = mul64s(s,s,p); }
    return r;
}

