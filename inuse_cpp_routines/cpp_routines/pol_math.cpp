#include<iostream> 
#include<cstdint> 
#include<random> 
#include<vector> 
#include<unordered_map>
#include<algorithm>
#include<time.h>
#include<chrono>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include"int128g.c"

using namespace std;
using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

ULNG seed=1;
ULNG mult=6364136223846793003LL;

/* 
INTEGER MATH ROUTINES: 
1. RAND64S
2. ADD64S 
3. SUB64S
4. MUL64S
5. MUL64ASM1
6. MUL64ASM2 
7. POWMOD64S
8. INVMOD64S

We are not using MUL64ASM1 and ASM2 routines. 
*/

/*
Function: rand64s
What it does:
  Generates one deterministic pseudo-random residue modulo p using the file's global 64-bit linear-congruential state.
Inputs:
  - p: modulus; the returned value is reduced modulo p.
Outputs:
  - Returns a residue in the range 0..p-1 and advances the global seed.
Example:
  LONG r = rand64s(101);  // one pseudo-random element of GF(101)
*/
LONG rand64s(LONG p){
    LONG x,y;
    extern ULNG seed,mult;
    seed=mult*seed;
    x=seed>>32;
    seed=mult*seed;
    y=seed>>32;
    x=(x<<31) | y;
    x=x%p;
    return(x);
}

/*
Function: add64b
What it does:
  Adds two finite-field residues and reduces the result modulo p without using division.
Inputs:
  - a, b: residues, normally assumed to satisfy 0 <= a,b < p.
  - p: modulus.
Outputs:
  - Returns (a+b) mod p.
Example:
  LONG c = add64b(70, 50, 101);  // c = 19
*/
inline LONG add64b(LONG a,LONG b,LONG p){
    LONG r=(a+b)-p;
    r+=(r>>63)&p;
    return r;
}

/*
Function: sub64b
What it does:
  Subtracts two finite-field residues and reduces the result modulo p.
Inputs:
  - a, b: residues, normally assumed to satisfy 0 <= a,b < p.
  - p: modulus.
Outputs:
  - Returns (a-b) mod p in the standard nonnegative residue range.
Example:
  LONG c = sub64b(3, 5, 101);  // c = 99
*/
inline LONG sub64b(LONG a,LONG b,LONG p){
    LONG r=(a-b);
    r+=(r>>63)&p;
    return r;
}

/*
Function: mul64b
What it does:
  Multiplies two residues with a 128-bit intermediate so the full product can be reduced safely modulo p.
Inputs:
  - a, b: residues to multiply.
  - p: modulus.
Outputs:
  - Returns (a*b) mod p.
Example:
  LONG c = mul64b(25, 9, 101);  // c = 23
*/
inline LONG mul64b(LONG a,LONG b, LONG p){
    ULNG128 res=(ULNG128)a*b;
    ULNG r=(ULNG)(res%p);
    return r;
}

/*
Function: neg64s
What it does:
  Computes the additive inverse of a residue modulo p.
Inputs:
  - a: residue in 0..p-1.
  - p: modulus.
Outputs:
  - Returns 0 when a=0; otherwise returns p-a.
Example:
  LONG b = neg64s(7, 101);  // b = 94
*/
inline LONG neg64s(LONG a,LONG p){ 
    return (a==0)?0:p-a; 
};

// We are assuming 0<=a,b<p for the following routines. 

/*
Function: mul64bASM
What it does:
  Computes modular multiplication using x86-64 mulq/divq instructions for the product and remainder.
Inputs:
  - a, b: residues to multiply.
  - p: modulus; inputs are assumed compatible with the assembly routine.
Outputs:
  - Returns (a*b) mod p.
Example:
  LONG c = mul64bASM(25, 9, 101);  // c = 23
*/
inline LONG mul64bASM(LONG a,LONG b,LONG p){
    LONG q, r;
    __asm__ __volatile__(           \
    "       mulq    %%rdx           \n\t" \
    "       divq    %4              \n\t" \
    : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
    return r;
}

/*
Function: mul64bASM2
What it does:
  Alternative x86-64 assembly implementation of modular multiplication using an explicit scratch register for p.
Inputs:
  - a, b: residues to multiply.
  - p: modulus.
Outputs:
  - Returns (a*b) mod p.
Example:
  LONG c = mul64bASM2(25, 9, 101);  // c = 23
*/
inline LONG mul64bASM2(LONG a,LONG b,LONG p){
    LONG q;
    LONG r;
    __asm__ __volatile__(
        "movq %[p],%%r8\n\t"
        "movq %[a],%%rax\n\t"
        "mulq %[b]\n\t"
        "divq %%r8\n\t"
        :"=&a"(q),"=&d"(r)
        :[a] "r"(a),[b] "rm"(b),[p] "r"(p)
        :"r8","cc"
    );
    return r;
}

/*
Function: powmod64s
What it does:
  Evaluates a^n modulo p using the modular multiplication routines and the exponentiation loop implemented here.
Inputs:
  - a: base; a negative representative is shifted into the field range.
  - n: nonnegative exponent.
  - p: modulus.
Outputs:
  - Returns the residue produced for a^n mod p; n=0 returns 1.
Example:
  LONG r = powmod64s(3, 5, 101);
*/
inline LONG powmod64s(LONG a,LONG n,LONG p){   
    LONG r,s;
    a+=(a>>63)&p;
    if(n==0){return 1;}
    if(n==1){return a;}
    for(r=1,s=a;n>0;n/=2){ 
        if(n&1){
            r=mul64b(r,s,p); 
            s=mul64b(s,s,p); 
        }
    }
    return r;
};

/*
Function: modinv64b
What it does:
  Computes the multiplicative inverse of c modulo p with the extended Euclidean algorithm.
Inputs:
  - c: value whose inverse is requested.
  - p: modulus.
Outputs:
  - Returns c^(-1) mod p when gcd(c,p)=1; returns 0 when no inverse exists.
Example:
  LONG inv = modinv64b(3, 101);  // inv = 34
*/
inline LONG modinv64b(LONG c,LONG p){   
    LONG d,r,q,r1,c1,d1;
    d=p;
    c1=1;
    d1=0;
    while(d!=0){
        q=c/d;
        r=c-q*d; 
        r1=c1-q*d1;
        c=d;
        c1=d1;
        d=r; 
        d1=r1;
    }
    if(c!=1) return(0);
    if(c1<0) c1+=p;
    return c1;
};

/*
POLYNOMIAL ARITHMETHIC ROUTINES: 
*/

struct RatReconFastWS{
    vector<LONG> r1;
    vector<LONG> r2;
    vector<LONG> t1;
    vector<LONG> t2;
    vector<LONG> q;
    vector<LONG> tmpT;

    /*
    Function: RatReconFastWS
    What it does:
      Constructs and sizes the reusable workspace used by fast rational reconstruction, avoiding repeated allocations inside the Euclidean loop.
    Inputs:
      - degM: maximum degree used to size all remainder, cofactor, quotient, and temporary arrays.
    Outputs:
      - Creates a RatReconFastWS object whose vectors have degM+1 slots; t1, t2, q, and tmpT are zero-initialized.
    Example:
      RatReconFastWS W(20);  // workspace for arrays with up to 21 coefficients
    */
    RatReconFastWS(int degM){
        int n=degM+1;
        r1.resize(n);
        r2.resize(n);
        t1.assign(n,0);
        t2.assign(n,0);
        q.assign(n,0);
        tmpT.assign(n,0);
    }
};

struct GCDEX{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;	
	int degS;
	int degT;
};

struct pairRFR{
	vector<LONG> r;
	vector<LONG> t;
	int degR;
	int degT;
	int flag;
};

struct GCDEXHIST{
	GCDEX g;
	vector<vector<LONG>> rTrace;
	vector<vector<LONG>> sTrace;
	vector<vector<LONG>> tTrace;
	vector<int> degRT;
	vector<int> degST;
	vector<int> degTT;
};

// Fast CPU routines 

#define ZMUL(z,a,b) do { \
    __asm__( \
    "       mulq    %%rdx   \n\t" \
            : "=a"(z[0]), "=d"(z[1]) \
            : "a"(a), "d"(b)); \
} while (0)

#define ZFMA(z,a,b) do { \
    unsigned long u,v; \
    __asm__( \
    "       mulq    %%rdx           \n\t" \
    "       addq    %%rax, %0       \n\t" \
    "       adcq    %%rdx, %1       \n\t" \
            : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) \
            : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b)); \
} while (0)

#define ZMOD(z,p) __asm__(\
    "       divq    %4      \n\t" \
    "       xorq    %0, %0  \n\t" \
            : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1] < p ? z[1] : z[1] % p), "r"(p))

/*
Function: genVEC64
What it does:
  Builds a random polynomial coefficient vector over GF(p) with coefficients generated by rand64s.
Inputs:
  - deg: requested polynomial degree; -1 represents the zero polynomial.
  - p: modulus used for each random coefficient.
Outputs:
  - Returns an empty vector for deg=-1, otherwise a vector of deg+1 residues.
Example:
  auto a = genVEC64(3, 101);  // four coefficients a[0]..a[3]
*/
vector<LONG> genVEC64(const int deg,const LONG p){
vector<LONG> v;
/*
Time complexity: O(d+1). 
Space complexity: O(d+1).
Auxillary space: O(1).
*/
if(deg==-1){return v;}
v.resize(deg+1);
for(int i=0;i<=deg;i++){
     v[i]=rand64s(p);
}
return v;
};

/*
Function: vecCOPY64
What it does:
  Creates an independent copy of a coefficient vector.
Inputs:
  - v: source vector.
Outputs:
  - Returns a new vector containing the same elements as v.
Example:
  auto b = vecCOPY64(a);
*/
vector<LONG> vecCOPY64(const vector<LONG> &v){
vector<LONG> temp; 
temp=v;
return temp;
};

/*
Function: vecfill64s
What it does:
  Fills a raw array with one constant value.
Inputs:
  - x: value to write.
  - A: destination array with at least n elements.
  - n: number of entries to fill.
Outputs:
  - Returns nothing; overwrites A[0]..A[n-1] with x.
Example:
  LONG A[4]; vecfill64s(0, A, 4);
*/
void vecfill64s( LONG x, LONG *A, int n )
{   int i;
    for( i=0; i<n; i++ ) A[i] = x;
    return;
}

/*
Function: polcopy64s
What it does:
  Copies the coefficients of a polynomial represented in ascending degree order.
Inputs:
  - A: source coefficients A[0]..A[d].
  - d: degree of the source polynomial.
  - B: destination array with at least d+1 slots.
Outputs:
  - Returns nothing; writes a copy of A into B.
Example:
  LONG A[3]={1,2,3}, B[3]; polcopy64s(A,2,B);
*/
void polcopy64s( LONG *A, int d, LONG *B )
{   int i;
    for( i=0; i<=d; i++) B[i]=A[i];
    return;
}

/*
Function: monic64s
What it does:
  Normalizes a nonzero polynomial over GF(p) so its leading coefficient becomes 1.
Inputs:
  - A: coefficient array, modified in place.
  - d: degree of A; negative means the zero polynomial.
  - p: modulus.
Outputs:
  - Returns nothing; when d>=0, A is scaled by the inverse of A[d] unless already monic.
Example:
  LONG A[3]={2,4,3}; monic64s(A,2,101);
*/
void monic64s(LONG *A,int d,LONG p) {
    int i; LONG inv;
    if(d<0 || A[d]==1) return;
    inv = modinv64b(A[d],p);
    for(i=0; i<d; i++) A[i] = mul64b(inv,A[i],p);
    A[d] = 1;
    return;
}

/*
Function: dispVEC64
What it does:
  Prints a coefficient vector as a human-readable polynomial in x, using ascending coefficient storage.
Inputs:
  - v: polynomial coefficient vector.
Outputs:
  - Returns nothing; writes the formatted polynomial to standard output.
Example:
  dispVEC64(vector<LONG>{1,2,3});  // prints 1*x^0 + 2*x^1 + 3*x^2
*/
void dispVEC64(const vector<LONG> &v){
if(v.size()==0) cout<<"O"<<"\n";
cout<<"[ ";
for(int i=0;i<v.size();i++){
    if(i==v.size()-1){
        cout<<v[i]<<"*x^"<<i<<"";
        break;
    }
    cout<<v[i]<<"*x^"<<i<<" + "<<"";
    }
    cout<<" ]"<<'\n';
};

// Returns a pair containing the new vector c=(a+b) mod p
// and the degree of c. 

/*
Function: pADDNEW64
What it does:
  Adds two polynomials over GF(p) and stores the sum in a newly allocated vector.
Inputs:
  - a, b: input coefficient vectors.
  - degA, degB: their degrees; -1 denotes the zero polynomial.
  - p: modulus.
Outputs:
  - Returns {c,degC}, where c=(a+b) mod p is trimmed to its true degree.
Example:
  auto [c,dc] = pADDNEW64(a,b,da,db,101);
*/
pair<vector<LONG>,int> pADDNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p){
	vector<LONG> c;
	int degC=-1;
	if(degA==-1 && degB==-1) return{c,degC};
	if(degA==-1)return{b,degB};
	if(degB==-1)return{a,degA};
	degC=max(degA,degB);
	c.resize(degC+1,0);
	int i=0;
	while(i<=degA && i<=degB){
		c[i]=add64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degA;i++){
		c[i]=a[i];
	}
	for(;i<=degB;i++){
		c[i]=b[i];
	}
	while(degC>=0 && c[degC]==0) degC--;
	if(degC==-1){
		c.clear();
		return{c,degC};
	}
	c.resize(degC+1);
	return{c,degC};
}

// In place addition. Overwrites a and returns the new degree.

/*
Function: pADDIP64
What it does:
  Adds polynomial b to polynomial a modulo p, modifying a rather than allocating a separate result.
Inputs:
  - a: first polynomial and destination.
  - b: second polynomial.
  - degA: degree of a, passed by reference and updated.
  - degB: degree of b.
  - p: modulus.
Outputs:
  - Returns the new degree of a and updates both a and degA.
Example:
  int dc = pADDIP64(a,b,da,db,101);
*/
int pADDIP64(vector<LONG> &a,const vector<LONG> &b,int &degA,const int degB,const LONG p){
	if(degA==-1&&degB==-1) return -1;
	if(degB==-1) return degA;
	if(degA==-1){
		a=b;
		degA=degB;
		return degA;
	}
	int maxDeg=max(degA,degB);
	if(a.size()<maxDeg+1) a.resize(maxDeg+1,0);
	int i=0;
	while(i<=degA && i<=degB){
		a[i]=add64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degB;i++) a[i]=b[i];
	while(maxDeg>=0 && a[maxDeg]==0) maxDeg--;
	if(maxDeg==-1){
		a.clear();
		return maxDeg;
	}
	a.resize(maxDeg+1);
	degA=maxDeg;
	return degA;
}

// Returns a pair containing the new vector c=(a-b) mod p
// and the degree of c. 

/*
Function: pSUBNEW64
What it does:
  Computes a-b over GF(p) and stores the result in a new vector.
Inputs:
  - a, b: input coefficient vectors.
  - degA, degB: their degrees; -1 denotes zero.
  - p: modulus.
Outputs:
  - Returns {c,degC}, where c=(a-b) mod p is trimmed to its true degree.
Example:
  auto [c,dc] = pSUBNEW64(a,b,da,db,101);
*/
pair<vector<LONG>,int> pSUBNEW64(const vector<LONG> &a,const vector<LONG> &b,
                                const int degA,const int degB,const LONG p){
    vector<LONG> c;
    int degC=-1;
    if(degA==-1 && degB==-1) return{c,degC};
    if(degA==-1){
        c.resize(degB+1,0);
        for(int i=0;i<=degB;i++) c[i]=neg64s(b[i],p);
        return{c,degB};
    }
    if(degB==-1) return{a,degA};

    degC=max(degA,degB);
    c.resize(degC+1,0);

    int i=0;
    while(i<=degA && i<=degB){
        c[i]=sub64b(a[i],b[i],p);
        i++;
    }
    for(; i<=degA; i++) c[i]=a[i];
    for(; i<=degB; i++) c[i]=neg64s(b[i],p);

    while(degC>=0 && c[degC]==0) degC--;
    if(degC==-1){ c.clear(); return{c,degC}; }
    c.resize(degC+1);
    return{c,degC};
}


// In place subtraction. Overwrites a and returns the new degree.

/*
Function: pSUBIP64
What it does:
  Subtracts polynomial b from the raw coefficient array a in place.
Inputs:
  - a: destination array with enough storage for the larger degree.
  - b: polynomial to subtract.
  - degA, degB: input degrees.
  - p: modulus.
Outputs:
  - Returns the degree of the updated a; -1 means the result is zero.
Example:
  int dc = pSUBIP64(a,b,da,db,101);
*/
int pSUBIP64(LONG *a,
             const LONG *b,
             int degA,
             const int degB,
             const LONG p){
	if(degA==-1&&degB==-1) return -1;
	if(degB==-1) return degA;
    if(degA==-1){
		for(int i=0;i<=degB;i++){
			a[i]=neg64s(b[i],p);
		}
		return degB;
	}
	int maxDeg=max(degA,degB);
	int i=0;
	while(i<=degA&&i<=degB){
		a[i]=sub64b(a[i],b[i],p);
		i++;
	}
	for(;i<=degB;i++) a[i]=neg64s(b[i],p);
	while(maxDeg>=0&&a[maxDeg]==0) maxDeg--;
	
	return maxDeg;
}

/*
Function: polsub64s
What it does:
  Computes the polynomial difference A-B over GF(p) into a separate output array C.
Inputs:
  - a, b: source coefficient arrays.
  - c: destination array.
  - da, db: source degrees.
  - p: modulus.
Outputs:
  - Returns degree(C) after trimming leading zeros.
Example:
  int dc = polsub64s(A,B,C,da,db,101);
*/
int polsub64s(LONG *a, LONG *b, LONG *c, int da, int db, LONG p) {
       int i,m;
       m = min(da,db);
       for( i=0; i<=m; i++ ) c[i] = sub64b(a[i],b[i],p);
       if( da==db ) { while( da>=0 && c[da]==0 ) da--; return da; }
       if( da>db ) { for ( i=db+1; i<=da; i++ ) c[i] = a[i]; return da; }
       for( i=da+1; i<=db; i++ ) c[i] = neg64s(b[i],p); return db;
    }

// Returns a pair containing the new vector c=(a*b) mod p
// and the degree of c. 

/*
Function: pMULNEW64
What it does:
  Multiplies two polynomials over GF(p) by coefficient convolution and allocates a new result vector.
Inputs:
  - a, b: input coefficient vectors.
  - degA, degB: their degrees.
  - p: modulus.
Outputs:
  - Returns {c,degC}; a zero input gives an empty vector with degree -1.
Example:
  auto [c,dc] = pMULNEW64(a,b,da,db,101);
*/
pair<vector<LONG>,int> pMULNEW64(const vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p){
	vector<LONG> c;
	if(degA<0 || degB<0) return {c,-1};
	int degC=degA+degB;
	c.resize(degC+1,0);
	for(int i=0;i<=degA;i++){
		for(int j=0;j<=degB;j++){
			LONG prod=mul64bASM(a[i],b[j],p);
			c[i+j]=add64b(c[i+j],prod,p);
		}
	}
	while(degC>=0 && c[degC]==0) degC--;
	if(degC==-1){
		c.clear();
		return {c,degC};
	}
	return {c,degC};
}

// In place multiplication. Overwrites a and returns the new degree.

/*
Function: pMULIP64
What it does:
  Multiplies a by b over GF(p), writing the product back into a and using fast accumulator paths selected by the size of p.
Inputs:
  - a: first factor and destination; caller must provide degA+degB+1 slots.
  - b: second factor.
  - degA, degB: input degrees.
  - p: modulus.
Outputs:
  - Returns the true degree of the product stored in a.
Example:
  int dc = pMULIP64(A,B,da,db,101);
*/
int pMULIP64(LONG *a,
             const LONG* b,
             int degA,
             int degB,
             const LONG p){
	if(degA<0 || degB<0) return -1;
	int i;
	int k;
	int m;
	int degC=degA+degB; // Called must guarantee that a has enough storage.
	/* 
	If p<2^31 then our product fits inside 2^63 bits.
	We essentially perform a convolution i.e. sum over i of
	a[i]*b[k-i].
	*/
	if(p<2147483648LL){
		LONG t;
		LONG p2; 
		p2=p<<32;
		for(k=degC;k>=0;k--){
			i=max(0,k-degB);
			m=min(k,degA);
			t=0LL; // Accumalator for the running sum.
			while(i<m){
				t-=a[i]*b[k-i];
				i++; 
				t-=a[i]*b[k-i];
				i++;
				t+=(t>>63)&p2;
			}
			if(i==m) t-=a[i]*b[k-i];
			t=(-t)%p;
			t+=(t>>63)&p;
			a[k]=t;
		}
	}
	else{
		ULNG z[2];
		for(k=degC;k>=0;k--){
			i=max(0,k-degB);
			m=min(k,degA);
			z[0]=z[1]=0LL; // 128 bit accumalators.
			while(i<m){
				//
				ZFMA(z,a[i],b[k-i]);
				i++;
				ZFMA(z,a[i],b[k-i]);
				i++;
				if(z[1]>=p) z[1]-=p;
			}
			if(i==m) ZFMA(z,a[i],b[k-i]);
			ZMOD(z,p);
			a[k]=z[0];
		}
	}
	while(degC>=0 && a[degC]==0) degC--;
	return degC;
}

/*
Function: polMUL64P
What it does:
  Reference in-place polynomial multiplication routine that computes a <- a*b modulo p with straightforward modular products.
Inputs:
  - a: first factor and destination with enough storage.
  - b: second factor.
  - degA, degB: degrees.
  - p: modulus.
Outputs:
  - Returns the degree of the product, or -1 for a zero factor.
Example:
  int dc = polMUL64P(A,B,da,db,101);
*/
int polMUL64P(LONG *a,
              LONG *b,
              int degA,
              int degB,
              LONG p){
    int i;
    int k;
    int m;
    LONG t;
    if(degA<0 || degB<0){
        return -1;
    }
    int degC=degA+degB;
    for(k=degC;k>=0;k--){
        i=max(0,k-degB);
        m=min(k,degA);
        for(t=0;i<=m;i++){
            t=add64b(t,mul64b(a[i],b[k-i],p),p);
        }
        a[k]=t;
    }
    while(degC>=0 && a[degC]==0){
        degC--;
    }
    return degC;
}

/*
Function: polmul64s
What it does:
  Multiplies two raw polynomial arrays over GF(p) into a separate array, using 64-bit or 128-bit accumulation depending on p.
Inputs:
  - A, B: source coefficient arrays.
  - C: destination with da+db+1 slots.
  - da, db: source degrees.
  - p: modulus.
Outputs:
  - Returns the degree of C, or -1 if either input is zero.
Example:
  int dc = polmul64s(A,B,C,da,db,101);
*/
int polmul64s( LONG * A, LONG * B, LONG * C, int da, int db, LONG p)
{
    int i,k,m;
    if( da<0 || db<0 ) return -1;
    int dc = da+db;
if( p<2147483648ll ) { LONG t, p2;
    p2 = p<<32;
    for( k=dc; k>=0; k-- ) {
       i = max(0,k-db);
       m = min(k,da);
       t = 0ll;
       while( i<m ) {
           t -= A[i]*B[k-i]; i++;
           t -= A[i]*B[k-i]; i++;
           t += (t>>63) & p2;
       }
       if( i==m ) t -= A[i]*B[k-i];
       t = (-t) % p;
       t += (t>>63) & p;
       C[k] = t;
    }
} else {ULNG z[2];
    for( k=dc; k>=0; k-- ) {
       i = max(0,k-db);
       m = min(k,da);
       z[0] = z[1] = 0ll;
       while( i<m ) {
           ZFMA(z,A[i],B[k-i]); i++;
           ZFMA(z,A[i],B[k-i]); i++;
           if( z[1]>=p ) z[1] -= p;
       }
       if( i==m ) ZFMA(z,A[i],B[k-i]);
       ZMOD(z,p);
       C[k] = z[0];
    }
}
    for( ; dc>=0 && C[dc]==0; dc-- );
    return( dc );
}

/*
Function: polfms64s
What it does:
  Performs a fused polynomial multiply-subtract C <- C - A*B modulo p.
Inputs:
  - A, B: multiplicand coefficient arrays.
  - C: polynomial updated in place.
  - da, db, dc: current degrees of A, B, and C.
  - p: modulus.
Outputs:
  - Returns the new degree of C after the subtraction.
Example:
  dc = polfms64s(A,B,C,da,db,dc,101);
*/
int polfms64s(LONG *A, LONG *B, LONG *C, int da, int db, int dc, LONG p)
{   // polynomial fused multiply subtract: C -= A*B
    int i,k,m; ULNG z[2];
    if( da<0 || db<0 ) return dc;

    for( k=0; k<=da+db; k++ ) {
        i = max(0, k-db);
        m = min(k, da);
        z[0] = z[1] = 0ll;

        while( i<m ) {
            ZFMA(z, A[i],   B[k-i]); i++;
            ZFMA(z, A[i],   B[k-i]); i++;
            if( z[1] >= p ) z[1] -= p;
        }
        if( i==m ) ZFMA(z, A[i], B[k-i]);

        ZMOD(z, p);

        if( k > dc ) {
            C[k] = (z[0] == 0 ? 0 : p - z[0]);
        } else {
            C[k] = sub64b(C[k], z[0], p);
        }
    }

    for( dc=max(dc, da+db); dc>=0 && C[dc]==0; dc-- );
    return dc;
}

/*
Function: pMULIP64VANDER
What it does:
  Multiplies two vector-backed polynomials for the Vandermonde code, placing a*b into b while normalizing coefficients modulo p.
Inputs:
  - a: left factor.
  - b: right factor and destination.
  - degA, degB: input degrees.
  - p: modulus.
Outputs:
  - Returns degree(b) after multiplication; b is replaced by the product.
Example:
  int dc = pMULIP64VANDER(a,b,da,db,101);
*/
int pMULIP64VANDER(vector<LONG> &a, vector<LONG> &b, int degA, int degB, const LONG p) {
    // Computes: b <- a * b   (in-place on b)
    // a is treated as the left factor (degree degA)
    // b is treated as the right factor and output (degree degB -> degA+degB)

    if (degA < 0 || degB < 0) {
        b.clear();
        return -1;
    }

    if ((int)a.size() < degA + 1) a.resize(degA + 1, 0);
    if ((int)b.size() < degB + 1) b.resize(degB + 1, 0);

    /*
    Function: normp (local lambda)
    What it does:
      Normalizes one signed coefficient into the standard residue range 0..p-1 for pMULIP64VANDER.
    Inputs:
      - x: coefficient to normalize. The surrounding lambda captures p by value.
    Outputs:
      - Returns x mod p as a nonnegative residue.
    Example:
      LONG r = normp(-1);  // r = p-1
    */
    auto normp = [p](LONG x) -> LONG {
        x %= p;
        if (x < 0) x += p;
        return x;
    };

    // Copy inputs (safe even if someone accidentally aliases a and b)
    vector<LONG> A(a.begin(), a.begin() + degA + 1);
    vector<LONG> B(b.begin(), b.begin() + degB + 1);

    for (int i = 0; i <= degA; ++i) A[i] = normp(A[i]);
    for (int j = 0; j <= degB; ++j) B[j] = normp(B[j]);

    int degC = degA + degB;
    vector<LONG> C(degC + 1, 0);

    // Convolution: C[k] = sum_{i=0}^degA A[i] * B[k-i]
    for (int i = 0; i <= degA; ++i) {
        if (A[i] == 0) continue;
        for (int j = 0; j <= degB; ++j) {
            if (B[j] == 0) continue;
            C[i + j] = add64b(C[i + j], mul64bASM(A[i], B[j], p), p);
        }
    }

    // Trim trailing zeros
    while (degC >= 0 && C[degC] == 0) --degC;

    if (degC < 0) {
        b.assign(1, 0);
        return -1;
    }

    C.resize(degC + 1);
    b.swap(C);
    return degC;
}

// Computes scalar polynomial multiplication i.e.
// A=c*A(x) where c is some scalar and returns a new vector.

/*
Function: polSCMULNEW64
What it does:
  Multiplies a polynomial by a scalar modulo p and returns a new coefficient vector.
Inputs:
  - a: source polynomial.
  - x: scalar multiplier; the special values 1 and -1 are handled directly.
  - degA: polynomial degree.
  - p: modulus.
Outputs:
  - Returns the scaled polynomial without changing a.
Example:
  auto b = polSCMULNEW64(a,5,da,101);
*/
vector<LONG> polSCMULNEW64(vector<LONG> &a,LONG x,int degA,const LONG p){
	vector<LONG> temp;
	if(x==1){return a;}
	temp.resize(degA+1,0);
	if(x==-1){
		for(int i=0;i<=degA;i++){
			temp[i]=neg64s(a[i],p);
		}
		return temp;
	}
	else{
		for(int i=0;i<=degA;i++){
			temp[i]=mul64bASM(a[i],x,p);
		}
	}
	return temp;	
}

// Computes scalar polynomial multiplication in place i.e.
// A=c*A(x) where c is some scalar. 

/*
Function: polSCMULIP64
What it does:
  Scales a polynomial by a scalar modulo p in place.
Inputs:
  - a: polynomial modified in place.
  - x: scalar multiplier.
  - degA: polynomial degree.
  - p: modulus.
Outputs:
  - Returns nothing; a becomes x*a mod p.
Example:
  polSCMULIP64(a,5,da,101);
*/
void polSCMULIP64(vector<LONG> &a,LONG x,int degA,const LONG p){
	// Since scalar value is 1 no difference.
	if(x==1){return;}
	// If it is negative we just use the neg function.
	if(x==-1){
		for(int i=0;i<=degA;i++){
			a[i]=neg64s(a[i],p);
		}
	}
	else{
		for(int i=0;i<=degA;i++){
			a[i]=mul64bASM(a[i],x,p);
		}
	}
}

// Computes A=A-(ax+b)*B efficiently using accumalators.

/*
Function: polSUBMUL64
What it does:
  Efficiently performs A <- A - (aVal*x+bVal)B over GF(p), a common degree-one Euclidean-algorithm step.
Inputs:
  - a: polynomial A, modified in place and with storage through index degB+1.
  - b: polynomial B.
  - aVal, bVal: coefficients of the linear quotient aVal*x+bVal.
  - degA, degB: current degrees.
  - p: modulus.
Outputs:
  - Returns the new degree of A.
Example:
  degA = polSUBMUL64(A,B,q1,q0,degA,degB,101);
*/
int polSUBMUL64(LONG *a,
                const LONG *b,
                LONG aVal,
                LONG bVal,
                int degA,
                int degB,
                const LONG p){
	ULNG z[2];
	LONG t;
	int i;
	/*
	If b is the zero polynomial, then A-(ax+b)*B=A so we simply
	return the degree of A.
	*/
	if(degB<0){
        return degA;
    }
	z[0]=z[1]=0LL;
	/*
	If degA<=degB, then we pad A with zereos. The caller needs to 
    guarantee a[0..degB+1] exists.
	*/
	while(degA<=degB){
    	++degA;
    	a[degA]=0;
    }
	/*
	Constant term is special in the sense b*B does not 
	have any effect on the degrees so we can compute A=b*B directly.
	*/
	t=mul64bASM(bVal,b[0],p);
	a[0]=sub64b(a[0],t,p);
	/*
	Basic for loop using 128 bit accumalators to compute 
	A[i]-(aVal*x-bVal)*B[i].
	*/
	for(i=1;i<=degB;i++){
        z[0]=z[1]=0ULL;
		ZMUL(z,aVal,b[i-1]);
		ZFMA(z,bVal,b[i]);
		ZMOD(z,p);
		t=a[i]-(LONG)z[0];
		a[i]=t+((t>>63)&p);
	}
	/*
	Here, we are updating the new leading coefficient.
	*/
	t=mul64bASM(aVal,b[degB],p);
	a[degB+1]=sub64b(a[degB+1],t,p);
	while(degA>=0 && (a[degA]==0 || a[degA]==p)){
        degA--;
    }
	return degA;
}

/*
Function: polSUBMUL64P
What it does:
  Straightforward modular-product version of A <- A - (aVal*x+bVal)B, used as a comparison/reference path.
Inputs:
  - a: polynomial A modified in place.
  - b: polynomial B.
  - aVal, bVal: linear quotient coefficients.
  - degA, degB: degrees.
  - p: modulus.
Outputs:
  - Returns the new degree of A; the routine prints FAIL if the degree unexpectedly does not drop.
Example:
  degA = polSUBMUL64P(A,B,q1,q0,degA,degB,101);
*/
int polSUBMUL64P(LONG *a,
                 const LONG *b,
                 LONG aVal,
                 LONG bVal,
                 int degA,
                 int degB,
                 const LONG p){
    LONG s;
    LONG t;
    int i;
    int d;
    if(degB==-1){
        return degA;
    }
    d=degA;
    while(degA<=degB){
        ++degA;
        a[degA]=0;
    }
    t=mul64b(bVal,b[0],p);
    a[0]=sub64b(a[0],t,p);
    for(i=1;i<=degB;i++){
        t=mul64b(aVal,b[i-1],p);
        t=add64b(t,mul64b(bVal,b[i],p),p);
        a[i]=sub64b(a[i],t,p);
    }
    t=mul64b(aVal,b[degB],p);
    a[degB+1]=sub64b(a[degB+1],t,p);
    while(degA>=0 && (a[degA]==0 || a[degA]==p)){
        degA--;
    }
    if(degA==d){
        printf("FAIL");
    }
    return degA;
}

// Evaluates a polynomial using Horners rule.

/*
Function: evalHORN64
What it does:
  Evaluates a vector-backed polynomial at alpha modulo p using Horner's rule.
Inputs:
  - a: coefficients in ascending degree order.
  - alpha: evaluation point.
  - p: modulus.
Outputs:
  - Returns a(alpha) mod p.
Example:
  LONG y = evalHORN64(a,7,101);
*/
LONG evalHORN64(vector<LONG>& a,LONG alpha,LONG p){
    LONG r = 0LL;
	for (int k=a.size();k-->0;){
        r=add64b(mul64bASM(r,alpha,p),a[k],p);
    }
    return r;
}

/*
Function: pEVAL64
What it does:
  Evaluates a raw-array polynomial at x modulo p using Horner's rule.
Inputs:
  - a: coefficient array.
  - d: degree; -1 represents zero.
  - x: evaluation point.
  - p: modulus.
Outputs:
  - Returns a(x) mod p, or 0 when d=-1.
Example:
  LONG y = pEVAL64(A,da,7,101);
*/
LONG pEVAL64(LONG *a,int d,LONG x,const LONG p){
	int i;
	LONG r;
	if(d==-1){return 0;}
	for(r=a[d],i=d-1;i>=0;i--){
		r=add64b(a[i],mul64bASM(x,r,p),p);
	}
	return r;
}

// My implementation of division. Same idea as PDIVIP64 but
// I am also returning the degree of both the quotient 
// and remainder.

/*
Function: pDIVDEG
What it does:
  Divides polynomial a by b over GF(p) in place, storing the remainder in the low coefficients and the quotient in the high portion of a.
Inputs:
  - a: dividend and combined remainder/quotient workspace.
  - b: divisor.
  - degA, degB: input degrees.
  - p: modulus.
Outputs:
  - Returns {degQ,degR}. The contents of a are overwritten by the division layout.
Example:
  auto [dq,dr] = pDIVDEG(a,b,da,db,101);
*/
pair<int,int> pDIVDEG(vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p){
	/* 
	Changes A in place.
	Suppose degA>=degB>=0. Then, time complexity is O((A-B+1)*(B+1)).
	*/
	if(degB<0){
		cout<<"DIV by 0.\n";
		exit(1);
	}
	if(degA<degB)return {-1,degA};
	if(degB == 0){
    	LONG b0 = b[0] % p; if(b0 < 0) b0 += p;
    	if(b0 == 0){ cout<<"DIV by 0.\n"; exit(1); }
    	LONG inv0 = modinv64b(b0, p);
    	for(int i=0;i<=degA;i++) a[i] = mul64bASM(a[i], inv0, p);
    	return {degA, -1};
	}
	LONG LTB=b[degB];
	LONG invLTB=modinv64b(LTB,p);
	for(int i=degA;i>=degB;i--){
		int k=i-degB;
		LONG LR=a[i];
		if(LR==0){
			a[i]=0; 
			continue;
		}
		LONG prod1=mul64bASM(LR,invLTB,p);
		for(int j=0;j<=degB;j++){
			LONG prod2=mul64bASM(prod1,b[j],p);
			a[k+j]=sub64b(a[k+j],prod2,p);
		}
		a[i]=prod1;
	}
	int degQ=degA-degB;
    while(degQ>0 && a[degB+degQ]==0) --degQ;
	if(degQ==0 && a[degB]==0) degQ=-1;
	int degR=degB-1;
    while(degR>0 && a[degR]==0) --degR;
	if(degR==0 && a[0]==0) degR=-1;
    return {degQ,degR};
}

// We divide a(x) by b(x) and put the remainder in the bottom 
// half of a so a[0...deg(b)-1] and quotient in the top half 
// so a[degb...dega] and return the degree of the remainder.

/*
Function: polDIVIP64
What it does:
  Fast raw-array polynomial division over GF(p). The dividend array is reused to store remainder and quotient.
Inputs:
  - a: dividend/workspace; after the call a[0..degB-1] contains the remainder and a[degB..degA] contains the quotient.
  - b: nonzero divisor.
  - degA, degB: degrees.
  - p: modulus.
Outputs:
  - Returns the degree of the remainder; -1 means exact division.
Example:
  int dr = polDIVIP64(A,B,da,db,101);
*/
int polDIVIP64(LONG *a,
               const LONG *b,
               int degA,
               int degB,
               const LONG p){
    int dq;
	int dr;
	int k;
	int j;
	int m; 
	LONG t;
	LONG inv;
    if(degB<0){ 
		cout<<"DIV BY 0.\n"; 
		exit(1); 
	}
    if(degA<degB) return degA;
	/* if (degB == 0) {
    // Divide by constant b0 (must be invertible mod p)
    LONG b0 = b[0] % p; if (b0 < 0) b0 += p;
    if (b0 == 0) { cout << "DIV BY 0.\n"; exit(1); }

    LONG inv0 = modinv64b(b0, p);
    for (int i = 0; i <= degA; i++) a[i] = mul64bASM(a[i], inv0, p);

    // remainder is 0
    return -1;
}*/

	/*
	Special case: If we have degA=degB and we have a monic
	divisor. Since, degrees are the same, the leading term 
	of A is the quotient (call it T). So we do A-(t*B) but 
	the first term is implicitly cancelled so we run a loop from 
	0 to degA-1 and update A accordingly. The degree of the remainder 
	will be degB-1 and we can then chop off the trailing zeroes.
	*/
    if(degA==degB && b[degB]==1){
        t=a[degA];
        for(k=0;k<degA;k++){
            if(b[k]){
				a[k]=sub64b(a[k],mul64bASM(t,b[k],p),p);
			}
		}
		for(dr=degA-1;dr>=0 && a[dr]==0;dr--);
        return dr;
    }
    dq=degA-degB;
    dr=degB-1;
	/* 
	We are working in Zp[x] so inverses are guaranteed to 
	exist. If LC(B) is monic then the inverse is simply 1. Else, 
	we compute the inverse.
	*/
    if(b[degB]==1)inv = 1; 
	else inv=modinv64b(b[degB],p);
	if(p<2147483648LL){ 
	LONG p2;
    p2=p<<32;
	/*
	Same idea as multiplication. If we know p<2^31 we can use 
	accumaltors to make things faster. The indices are 
	extrmely confusing.
	*/
    for(k=degA;k>=0;k--){
        t=a[k];
        m=min(dr,k);
        j=max(0,k-dq);
        while(j<m){
            t-=b[j]*a[k-j+degB]; 
			j++;
            t-=b[j]*a[k-j+degB]; 
			j++;
            t+=(t>>63)&p2;
        }
        if(j==m){
            t-=b[j]*a[k-j+degB];
        }
        t=t%p;
        t+=(t>>63)&p;
        if(k>=degB && inv!=1){
            t=mul64bASM(t,inv,p);
        }
        a[k]=t;
    }
} else{
	ULNG z[2];
    for(k=degA;k>=0;k--){
        z[0]=z[1]=0LL;
        m=min(dr,k);
        j=max(0,k-dq);
        while(j<m){
            ZFMA(z,b[j],a[k-j+degB]); 
			j++;
            ZFMA(z,b[j],a[k-j+degB]); 
			j++;
            if(z[1]>=p)z[1]-=p;
        }
        if(j==m){
            ZFMA(z,b[j],a[k-j+degB]);
        }
        ZMOD(z,p);
        t=a[k]-z[0];
        t+=(t>>63)&p;
        if(k>=degB && inv!=1){
            t=mul64bASM(t,inv,p);
        }
        a[k]=t;
    }
}
    while(dr>=0 && a[dr]==0){
        dr--;
    }
    return dr;
}

/*
Function: polDIVP
What it does:
  Reference polynomial-division implementation that overwrites a with the same remainder/quotient layout as polDIVIP64.
Inputs:
  - a: dividend/workspace.
  - b: divisor.
  - degA, degB: degrees.
  - p: modulus.
Outputs:
  - Returns degree of the remainder, or -1 for exact division.
Example:
  int dr = polDIVP(A,B,da,db,101);
*/
int polDIVP(LONG *a,
            LONG *b,
            int degA,
            int degB,
            LONG p){
    int degQ;
    int degR;
    int k;
    int j;
    int m;
    LONG t;
    LONG inv; 
    if(degB<0){
        printf("Div by 0\n");
        return -1;
    }
    if(degA<degB){
        return degA;
    }
    degQ=degA-degB;
    degR=degB-1;
    if(b[degB]==1){
        inv=1;
    }
    else{
        inv=modinv64b(b[degB],p);
    }
    for(k=degA;k>=0;k--){
        t=a[k];
        m=min(degR,k);
        j=max(0,k-degQ);
        for(t=a[k];j<=m;j++){
            t=sub64b(t,mul64b(b[j],a[k-j+degB],p),p);
        }
        if(k>=degB && inv!=1){
            t=mul64b(t,inv,p);
        }
        a[k]=t;
    }
    while(degR>0 &&  a[degR]==0){
        degR--;
    }
    return degR;
}

// Makes a polynomial monic. We can do this as we are working 
// in Zp[x] and this is a field so inverses exist.

/*
Function: polMAKEMONIC64
What it does:
  Normalizes a vector-backed polynomial over GF(p) by scaling all coefficients so the leading coefficient is 1.
Inputs:
  - a: polynomial modified in place.
  - p: modulus.
Outputs:
  - Returns nothing; a is unchanged if it is zero or already monic.
Example:
  polMAKEMONIC64(a,101);
*/
void polMAKEMONIC64(vector<LONG> &a,const LONG p){
	int degA=a.size()-1;
	if(degA<0 || a[degA]==1) return;
	LONG invTerm;
	invTerm=modinv64b(a[degA],p);
	for(int i=0;i<degA;i++){
		a[i]=mul64bASM(invTerm,a[i],p);
	}
	a[degA]=1;
}

// My version of computing the gcd(a(x),b(x)). This returns 
// the updated vector a with the GCD and its degree.

/*
Function: polGCDNEW64
What it does:
  Computes the monic polynomial gcd of a and b over GF(p) with the Euclidean algorithm.
Inputs:
  - a, b: polynomial vectors; they are used as mutable Euclidean workspaces.
  - degA, degB: degrees.
  - p: modulus.
Outputs:
  - Returns {g,degG}, where g is the monic gcd.
Example:
  auto [g,dg] = polGCDNEW64(a,b,da,db,101);
*/
pair<vector<LONG>,int> polGCDNEW64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p){
	// Division dominates. 
	// Space is O((degA+1)+(degB+1)).
	if(degA==-1){
		polMAKEMONIC64(b,p);
		return {b,degB};
	}
	if(degB==-1){
		polMAKEMONIC64(a,p);
		return {a,degA};
	}
	if(degA<degB){
		swap(a,b);
		swap(degA,degB);
	}
	while(degB!=-1){
		pair<int,int> QR=pDIVDEG(a,b,degA,degB,p);
        int degR=QR.second;
		a.swap(b);       
        degA=degB;
        degB=degR;
    }
	polMAKEMONIC64(a,p);
	return {a,degA};
}

/*
Function: ratReconNormal
What it does:
  Runs the extended-Euclidean rational-reconstruction kernel on M and U, stopping at the requested numerator degree and normalizing the denominator to be monic.
Inputs:
  - m, u: coefficient vectors for M and U.
  - degM, degU: their degrees.
  - N, D: requested numerator/denominator degree bounds used by the reconstruction interface.
  - p: modulus.
  - W: reusable RatReconFastWS workspace.
  - rOut/tOut and degROut/degTOut: output buffers and returned degrees.
Outputs:
  - Returns 0 on reconstruction success; negative codes indicate that the target row was not found or normalization failed. On success rOut/tOut contain numerator/denominator coefficients.
Example:
  int rc = ratReconNormal(M,U,degM,degU,N,D,p,W,num,&dn,den,&dd);
*/
int ratReconNormal(const vector<LONG> &m,
    const vector<LONG> &u,
    int degM,
    int degU,
    int N,
    int D,
    const LONG p,
    RatReconFastWS &W,
    LONG *rOut,
    int &degROut,
    LONG *tOut,
    int &degTOut){

// Copy inputs into workspace
std::copy_n(m.data(), degM + 1, W.r1.data());
std::copy_n(u.data(), degU + 1, W.r2.data());

// Initialize t-sequence:
// r1 = m, r2 = u
// t1 = 0, t2 = 1
W.t2[0] = 1;

int degA  = degM;
int degB  = degU;
int degT1 = -1;
int degT2 = 0;

// Ensure degA >= degB initially
if(degA < degB){
std::swap(W.r1, W.r2);
std::swap(degA, degB);
std::swap(W.t1, W.t2);
std::swap(degT1, degT2);
}

while(degB != -1){

// Stop at the first index k such that deg(r_k) == N
if(degB == N){
    degROut = degB;
    degTOut = degT2;

    std::copy_n(W.r2.data(), degROut + 1, rOut);
    std::copy_n(W.t2.data(), degTOut + 1, tOut);

    // Normalize so denominator is monic
    if(degTOut >= 0){
        LONG lc = tOut[degTOut];
        if(lc == 0){
            degROut = -1;
            degTOut = -1;
            return -30; // unexpected bad denominator
        }

        if(lc != 1){
            LONG lcInv = modinv64b(lc, p);
            if(lcInv == 0){
                degROut = -1;
                degTOut = -1;
                return -31; // inverse does not exist
            }

            for(int i = 0; i <= degROut; i++){
                rOut[i] = mul64bASM(rOut[i], lcInv,p);
            }
            for(int i = 0; i <= degTOut; i++){
                tOut[i] = mul64bASM(tOut[i], lcInv,p);
            }
        }
    }

    return 0;
}

LONG uInv, aVal, bVal;
int degR, degQ, degT;

// Special degree-1 quotient step
if(degB > 0 && degA - degB == 1){
uInv = modinv64b(W.r2[degB], p);

aVal = mul64bASM(W.r1[degA], uInv, p);

bVal = mul64bASM(aVal, W.r2[degB - 1], p);
bVal = mul64bASM(uInv, sub64b(W.r1[degA - 1], bVal, p), p);

degR = polSUBMUL64(W.r1.data(), W.r2.data(),
           aVal, bVal, degA, degB, p);

degT = polSUBMUL64(W.t1.data(), W.t2.data(),
           aVal, bVal, degT1, degT2, p);
}
else{
// Divide r1 by r2:
// quotient goes into high part of W.r1, remainder stays in low part
degR = polDIVIP64(W.r1.data(), W.r2.data(), degA, degB, p);
degQ = degA - degB;

for(int i = 0; i <= degQ; i++){
W.q[i] = W.r1[degB + i];
}

if(degT2 >= 0){
for(int i = 0; i <= degT2; i++){
W.tmpT[i] = W.t2[i];
}

int degTmpT = degT2;
degTmpT = pMULIP64(W.tmpT.data(), W.q.data(), degTmpT, degQ, p);
degT = pSUBIP64(W.t1.data(), W.tmpT.data(), degT1, degTmpT, p);
}
else{
degT = degT1;
}
}

if(degR < 0){
break;
}

// Shift:
// (r1,r2) <- (r2,r)
// (t1,t2) <- (t2,t)
std::swap(W.r1, W.r2);
degA = degB;
degB = degR;

std::swap(W.t1, W.t2);
int oldDegT2 = degT2;
degT2 = degT;
degT1 = oldDegT2;
}

degROut = -1;
degTOut = -1;
return -20;
};

/* HFTRFR and DFTRFR */

/*
Function: polGCDvec
What it does:
  Internal vector-based monic gcd helper used by the fault-tolerant rational reconstruction routines.
Inputs:
  - a, b: polynomial copies passed by value.
  - degA, degB: degrees.
  - p: modulus.
  - g: output vector.
Outputs:
  - Returns degree(g) and writes the monic gcd into g.
Example:
  int dg = polGCDvec(A,da,B,db,101,g);
*/
static int polGCDvec(std::vector<LONG> a,int degA,std::vector<LONG> b,int degB,LONG p,std::vector<LONG>&g){
    if(degA<0 && degB<0){ g.assign(1,1); return 0; }             // gcd(0,0) := 1
    if(degB<0){ g=a; g.resize(degA+1); monic64s(g.data(),degA,p); return degA; }
    if(degA<0){ g=b; g.resize(degB+1); monic64s(g.data(),degB,p); return degB; }

    int m=std::max(degA,degB)+2;              // polSUBMUL64 may write index degB+1
    a.resize(m,0); b.resize(m,0);
    LONG *C=a.data(), *D=b.data(), *R;
    int dr;
    if(degA<degB){ R=C; C=D; D=R; dr=degA; degA=degB; degB=dr; }
    while(1){
        if(degB>0 && degA-degB==1){                              // quotient q1*x + q0
            LONG u  = modinv64b(D[degB],p);
            LONG q1 = mul64b(C[degA],u,p);
            LONG q0 = mul64b(q1,D[degB-1],p);
            q0 = mul64b(u,sub64b(C[degA-1],q0,p),p);
            dr = polSUBMUL64(C,D,q1,q0,degA,degB,p);             // C -= (q1 x + q0) D
        } else {
            dr = polDIVIP64(C,D,degA,degB,p);                    // C low part := C mod D
        }
        if(dr<0){                                                // D divides C, so gcd = D
            g.assign(D,D+degB+1);
            monic64s(g.data(),degB,p);
            return degB;
        }
        R=C; C=D; D=R; degA=degB; degB=dr;
    }
}

/*
Function: polExactQuo
What it does:
  Computes the quotient A/B in the case where B is expected to divide A exactly.
Inputs:
  - A, B: input polynomial vectors.
  - degA, degB: degrees.
  - p: modulus.
  - Q: output quotient vector.
Outputs:
  - Returns degree(Q) and writes the quotient to Q.
Example:
  int dq = polExactQuo(A,da,B,db,101,Q);
*/
static int polExactQuo(const std::vector<LONG>&A,int degA,const std::vector<LONG>&B,int degB,LONG p,std::vector<LONG>&Q){
    std::vector<LONG> a=A; a.resize(degA+1);
    polDIVIP64(a.data(),B.data(),degA,degB,p);
    int degQ=degA-degB;
    Q.assign(a.begin()+degB,a.begin()+degB+degQ+1);
    while(degQ>=0 && Q[degQ]==0) degQ--;
    Q.resize(degQ<0?1:degQ+1);
    return degQ;
}

/*
Function: eeaStep
What it does:
  Performs one extended-Euclidean step on consecutive remainder/cofactor pairs, using the optimized linear-quotient path when possible.
Inputs:
  - Rp,R and degRp,degR: previous/current remainders.
  - Tp,T and degTp,degT: matching cofactors.
  - p: modulus.
  - rem,drem,nT,dnT: output remainder and next cofactor with their degrees.
Outputs:
  - Returns nothing; writes the next Euclidean remainder and cofactor through the output references.
Example:
  eeaStep(Rp,dRp,R,dR,Tp,dTp,T,dT,101,rem,drem,nT,dnT);
*/
static void eeaStep(const std::vector<LONG>&Rp,int degRp,const std::vector<LONG>&R,int degR,
                    const std::vector<LONG>&Tp,int degTp,const std::vector<LONG>&T,int degT,LONG p,
                    std::vector<LONG>&rem,int &drem,std::vector<LONG>&nT,int &dnT){
    int degQ=degRp-degR;
    if(degR>0 && degQ==1){                                         // Q = aVal*x + bVal
        LONG uInv=modinv64b(R[degR],p);
        LONG aVal=mul64bASM(Rp[degRp],uInv,p);
        LONG bVal=mul64bASM(aVal,R[degR-1],p);
        bVal=mul64bASM(uInv,sub64b(Rp[degRp-1],bVal,p),p);
        rem=Rp; rem.resize(degRp+1,0);                             // needs index degR+1 == degRp
        drem=polSUBMUL64(rem.data(),R.data(),aVal,bVal,degRp,degR,p);
        nT=Tp;  nT.resize(std::max(degTp,degT+1)+1,0);             // cofactor grows to degT+1
        dnT=polSUBMUL64(nT.data(),T.data(),aVal,bVal,degTp,degT,p);
    } else {                                                       // general division
        std::vector<LONG> a=Rp; a.resize(degRp+1);
        drem=polDIVIP64(a.data(),R.data(),degRp,degR,p);
        std::vector<LONG> Q(a.begin()+degR,a.begin()+degR+degQ+1);
        if(drem>=0) rem.assign(a.begin(),a.begin()+drem+1); else rem.assign(1,0);
        std::vector<LONG> QT=T; QT.resize(degT+degQ+1,0);
        int dQT=pMULIP64(QT.data(),Q.data(),degT,degQ,p);          // QT := T*Q
        nT=Tp; nT.resize(std::max(degTp,dQT)+1,0);
        dnT=pSUBIP64(nT.data(),QT.data(),degTp,dQT,p);             // nT := Tp - T*Q
    }
}


/*
Function: dftrfr
What it does:
  Implements deterministic fault-tolerant rational function reconstruction from the modular congruence U mod M with numerator, denominator, and error bounds.
Inputs:
  - M,U and degM,degU: reconstruction polynomials.
  - N,D: numerator and denominator degree bounds.
  - E: allowed error count.
  - p: modulus.
  - numOut/denOut and degNum/degDen: output buffers and degrees.
Outputs:
  - Returns 1 on successful validated reconstruction and 0 when the data/bounds do not admit one.
Example:
  int ok = dftrfr(M,U,degM,degU,N,D,E,101,num,dn,den,dd);
*/
static int dftrfr(const std::vector<LONG>&M,const std::vector<LONG>&U,int degM,int degU,
                  int N,int D,int E,LONG p,
                  LONG *numOut,int &degNum,LONG *denOut,int &degDen){
    degNum=-1; degDen=-1;
    if(degM <= N+D+2*E) return 0;                          // not enough points
    if(degU < 0){ numOut[0]=0; degNum=0; denOut[0]=1; degDen=0; return 1; } // U==0 -> 0/1
    if(E>0){ std::vector<LONG> gg; int dg=polGCDvec(M,degM,U,degU,p,gg); if(dg>N+E) return 0; } // gcd test (cannot fire when E=0)

    std::vector<LONG> Rp=M,R=U,Tp(1,0),T(1,1);            // Rp=M,R=U ; Tp=0,T=1
    int degRp=degM,degR=degU,degTp=-1,degT=0;
    while(degR>=0 && degT<=D+E){
        std::vector<LONG> rem,nT; int drem,dnT;
        eeaStep(Rp,degRp,R,degR,Tp,degTp,T,degT,p,rem,drem,nT,dnT);
        Rp=R; degRp=degR; R=rem; degR=drem;
        Tp=T; degTp=degT; T=nT; degT=dnT;
    }
    std::vector<LONG> d; int degd=polGCDvec(Rp,degRp,Tp,degTp,p,d);
    std::vector<LONG> numV,denV;
    int degn =polExactQuo(Rp,degRp,d,degd,p,numV);
    int degdn=polExactQuo(Tp,degTp,d,degd,p,denV);
    if(degn>N || degdn>D) return 0;                        // degree validation
    LONG iv=modinv64b(denV[degdn],p);                      // make den monic
    for(int i=0;i<=degn;i++)  numOut[i]=mul64b(numV[i],iv,p);
    for(int i=0;i<=degdn;i++) denOut[i]=mul64b(denV[i],iv,p);
    degNum=degn; degDen=degdn;
    return 1;
}

/*
Function: hftrfr
What it does:
  Runs the heuristic/gap phase of fault-tolerant rational reconstruction, extracting the candidate numerator/denominator row, its common bad factor Lambda, and the largest quotient-degree gap.
Inputs:
  - M,U and degM,degU: input polynomials.
  - p: modulus.
  - fOut,gOut,lamOut: output coefficient arrays.
  - degF,degG,degLam,qmax: output degrees and maximum quotient degree.
Outputs:
  - Returns 0 after filling the candidate data; gOut is normalized to be monic and Lambda is monic.
Example:
  hftrfr(M,U,degM,degU,101,f,dF,g,dG,lam,dL,qmax);
*/
static int hftrfr(const std::vector<LONG>&M,const std::vector<LONG>&U,int degM,int degU,LONG p,
                  LONG *fOut,int &degF,LONG *gOut,int &degG,LONG *lamOut,int &degLam,int &qmax){
    std::vector<LONG> rp=M,rc=U,tp(1,0),tc(1,1);
    int degrp=degM,degrc=degU,degtp=-1,degtc=0;
    std::vector<LONG> f=M,g(1,1); int degf=degM,degg=0;
    qmax=0;
    while(degrc>=0){
        int degQ=degrp-degrc;
        if(degQ>qmax){ qmax=degQ; f=rc; degf=degrc; g=tc; degg=degtc; }   // save gap row
        std::vector<LONG> rem,nT; int drem,dnT;
        eeaStep(rp,degrp,rc,degrc,tp,degtp,tc,degtc,p,rem,drem,nT,dnT);
        rp=rc; degrp=degrc; rc=rem; degrc=drem;
        tp=tc; degtp=degtc; tc=nT; degtc=dnT;
    }
    std::vector<LONG> Lam; int degLambda=polGCDvec(f,degf,g,degg,p,Lam);
    std::vector<LONG> fc,gc;
    int degfc=polExactQuo(f,degf,Lam,degLambda,p,fc);
    int deggc=polExactQuo(g,degg,Lam,degLambda,p,gc);
    LONG iv=modinv64b(gc[deggc],p);                        // make g monic
    for(int i=0;i<=degfc;i++)     fOut[i]=mul64b(fc[i],iv,p);
    for(int i=0;i<=deggc;i++)     gOut[i]=mul64b(gc[i],iv,p);
    for(int i=0;i<=degLambda;i++) lamOut[i]=Lam[i];        // Lam already monic
    degF=degfc; degG=deggc; degLam=degLambda;
    return 0;
}

/* CPP DFTRFR */
/*
Function: cppDFTRFR
What it does:
  C-compatible wrapper around dftrfr for Maple's external-function interface, including argument checks and output initialization.
Inputs:
  - M/U arrays with lengths and degrees.
  - N,D,E: reconstruction bounds.
  - p: modulus.
  - nOut,dOut: output buffers with their allocated lengths.
Outputs:
  - Returns 0 on success, 1 when reconstruction fails, and a negative value for invalid arguments.
Example:
  int rc = cppDFTRFR(3,2,M,2,1,U,1,1,0,101,2,num,2,den);
*/
extern "C" int cppDFTRFR(int mLen,int degM,const LONG *M,
                         int uLen,int degU,const LONG *U,
                         int N,int D,int E,const LONG p,
                         int nOutLen,LONG *nOut,
                         int dOutLen,LONG *dOut){
    if(!M||!U||!nOut||!dOut)                          return -1;
    if(degM<0||degU<-1)                               return -1;   // degU=-1 => U==0
    if(mLen<=0||uLen<=0||nOutLen<=0||dOutLen<=0)      return -1;
    if(nOutLen<N+1||dOutLen<D+1)                      return -1;
    for(int i=0;i<nOutLen;i++) nOut[i]=0;
    for(int i=0;i<dOutLen;i++) dOut[i]=0;
    std::vector<LONG> m(M,M+degM+1), u;
    if(degU>=0) u.assign(U,U+degU+1); else u.assign(1,0);
    int degNum=-1,degDen=-1;
    int st=dftrfr(m,u,degM,degU,N,D,E,p,nOut,degNum,dOut,degDen);
    if(st==0) return 1;                                    // reconstruction FAILED
    return 0;
}

/* CPP HFTRFR */
/*
Function: cppHFTRFR
What it does:
  C-compatible wrapper around hftrfr that exposes the candidate factors and qmax to Maple.
Inputs:
  - M/U arrays with lengths and degrees.
  - p: modulus.
  - fOut,gOut,lamOut: output buffers and lengths.
  - info: output metadata array; info[0] receives qmax.
Outputs:
  - Returns 0 on success and a negative value for invalid arguments.
Example:
  int rc = cppHFTRFR(mLen,dM,M,uLen,dU,U,p,fLen,f,gLen,g,lLen,lam,1,info);
*/
extern "C" int cppHFTRFR(int mLen,int degM,const LONG *M,
                         int uLen,int degU,const LONG *U,
                         const LONG p,
                         int fOutLen,LONG *fOut,
                         int gOutLen,LONG *gOut,
                         int lamOutLen,LONG *lamOut,
                         int infoLen,LONG *info){
    if(!M||!U||!fOut||!gOut||!lamOut||!info) return -1;
    if(degM<0||degU<0)                       return -1;
    if(mLen<=0||uLen<=0||fOutLen<degM+1||gOutLen<degM+1||lamOutLen<degM+1||infoLen<1) return -1;
    for(int i=0;i<fOutLen;i++)  fOut[i]=0;
    for(int i=0;i<gOutLen;i++)  gOut[i]=0;
    for(int i=0;i<lamOutLen;i++)lamOut[i]=0;
    for(int i=0;i<infoLen;i++)  info[i]=0;
    std::vector<LONG> m(M,M+degM+1), u(U,U+degU+1);
    int degF=-1,degG=-1,degLam=-1,qmax=0;
    hftrfr(m,u,degM,degU,p,fOut,degF,gOut,degG,lamOut,degLam,qmax);
    info[0]=qmax;
    return 0;
}

/* VANDERMONDE SOLVER ROUTINES: */

/*
Function: VandermondeSolve64s
What it does:
  Solves the transposed Vandermonde system used by sparse interpolation, constructing the master polynomial product_i(x-m[i]) as workspace.
Inputs:
  - m: distinct interpolation nodes, modified as workspace-compatible input.
  - y: right-hand-side values.
  - n: system size.
  - a: output coefficients.
  - M: workspace with n+1 slots.
  - shift: optional monomial shift correction.
  - p: modulus.
Outputs:
  - Returns nothing; writes the solved coefficient vector to a[0..n-1].
Example:
  VandermondeSolve64s(nodes,values,n,coeffs,workspace,0,101);
*/
void VandermondeSolve64s(LONG *m,LONG *y,int n,LONG *a,LONG *M,int shift,LONG p)
{  
    int i,j;
    LONG u,s,A[2];
    M[0] = 1;
    A[1] = 1;
    for( i=0; i<n; i++ ){ 
         A[0] = neg64s(m[i],p);
         polmul64s(A,M,M,1,i,p);
    }
    for( j=0; j<n; j++ ){
         A[0] = neg64s(m[j],p);
         i = polDIVIP64(M,A,n,1,p); // M = M / (x-m[j])
         for( i=0; i<n; i++ ) M[i] = M[i+1]; // move quotient to front of M
         u = pEVAL64(M,n-1,m[j],p);
         if( u==0 ) { printf("Roots are not distinct!\n"); return; }
         u = modinv64b(u,p);
         // compute a[j] = M dot y
         for( s=0,i=0; i<n; i++ ) s = add64b(s,mul64b(M[i],y[i],p),p);
         s = mul64b(u,s,p);
         //s = mul64b(u,s,p);
         if( shift!=0 ) {
             u = modinv64b(m[j],p);
             u = powmod64s(u,shift,p);
             //s = mul64b(u,s,p);
             s = mul64b(u,s,p);
         }
         a[j] = s;
         polmul64s(A,M,M,1,n-1,p); // restore M = M x (x-m[j])
    }
    return;
}

/*
Function: cppVSolve
What it does:
  C-compatible wrapper for VandermondeSolve64s that protects Maple inputs by copying the node/value arrays before the destructive kernel runs.
Inputs:
  - mIn,yIn: node and value arrays with lengths mLen,yLen.
  - shift: shift correction.
  - p: modulus.
  - aOut/outLen: output buffer and capacity.
Outputs:
  - Returns 0 on success or a negative argument-validation code; writes n coefficients to aOut.
Example:
  int rc = cppVSolve(n,m,n,y,0,101,n,a);
*/
extern "C" int cppVSolve(int mLen,
    const LONG *mIn,
    int yLen,
    const LONG *yIn,
    int shift,
    const LONG p,
    int outLen,
    LONG *aOut)
{

// INITIAL CHECKS:

if (!mIn || !yIn || !aOut) {
    return -1;
}
if (mLen <= 0 || yLen <= 0 || outLen <= 0) {
    return -2;
}
if (mLen != yLen) {
    return -3;
}
if (outLen < mLen) {
    return -4;
}

const int n = mLen;

// INITIALIZING OUTPUT ARRAY:

for (int i = 0; i < outLen; ++i) {
    aOut[i] = 0;
}

// local working copies since kernel overwrites m, y, and M
vector<LONG> m(mIn, mIn + n);
vector<LONG> y(yIn, yIn + n);

// M holds the master polynomial prod_i (x - m[i]) of degree n,
// so it needs n+1 coefficient slots.
vector<LONG> M(n + 1, 0);

VandermondeSolve64s(m.data(), y.data(), n, aOut, M.data(), shift, p);

return 0;
}

/*
BERLEKAMP MASSEY ROUTINES:
*/

/*
Function: BerlekampMassey64s
What it does:
  Computes a Berlekamp-Massey connection polynomial for a finite-field sequence using a half extended-Euclidean algorithm.
Inputs:
  - a: sequence of N residues.
  - N: sequence length.
  - L: output coefficient array.
  - W: scratch workspace.
  - p: modulus.
Outputs:
  - Returns degree(L); returns -1 when no nontrivial connection polynomial is produced. L is made monic.
Example:
  int d = BerlekampMassey64s(seq,N,L,W,101);
*/
int BerlekampMassey64s(LONG *a,int N,LONG *L,LONG *W,LONG p)
{
    // Input sequence a = [a1,a2,a3,...,aN]
    // Output polynomial Lambda(x) is written to L
    // Uses the half extended Euclidean algorithm
    int i,m,n,dr,dq,dr0,dr1,dv0,dv1,dt;
    LONG *r,*q,*r0,*r1,*v0,*v1,*t,u,A,b;
    while( N>0 && a[N-1]==0 ) N--; // ignore leading zeroes
    n = N/2;
    N = 2*n;
    if( N==0 ) return -1;
    m = N-1;
    // W is space for r0 = x^N and r1 of degree m and v0 and v1 of degree at most n
    r0 = W; r1 = r0+N+1; v0 = r1+N; v1 = v0+n+1;
    vecfill64s(0,r0,N); r0[N] = 1; dr0 = N;             // r0 = x^(2*n)
    for(i=0; i<N; i++) r1[i] = a[m-i];
    for(dr1=m; dr1>=0 && r1[dr1]==0; dr1--);            // r1 = sum(a[m-i]*x^i,i=0..m)
    if( dr1==-1 ) return -1;
    dv0 = -1;                                           // v0 = 0
    v1[0] = 1; dv1 = 0;                                 // v1 = 1
    while( n <= dr1 ) {
        if( dr1>0 && dr0-dr1==1 ) { // normal case
            u = modinv64b(r1[dr1],p);
            A = mul64b(r0[dr0],u,p);
            b = mul64b(A,r1[dr1-1],p);
            b = mul64b(u,sub64b(r0[dr0-1],b,p),p);             // quotient q = A x + b
            //dr = polsubmulP(r0,r1,A,b,dr0,dr1,p,P);            // r0 = r0 - (A x + b) r1
            dr = polSUBMUL64(r0,r1,A,b,dr0,dr1,p);            // r0 = r0 - (A x + b) r1
            //dt = polsubmulP(v0,v1,A,b,dv0,dv1,p,P);            // v0 = v0 - (A x + b) v1
            dt = polSUBMUL64(v0,v1,A,b,dv0,dv1,p);            // v0 = v0 - (A x + b) v1
        } else {
           dr = polDIVIP64(r0,r1,dr0,dr1,p);
           q = r0+dr1; dq = dr0-dr1;                           // q = quo(r0,r1)
           dt = polmul64s(q,v1,L,dq,dv1,p);
           dt = polsub64s(v0,L,v0,dv0,dt,p);
        }
        r = r0; r0 = r1; r1 = r; dr0 = dr1; dr1 = dr;         // r0,r1 = r1,rem(r0,r1)
        t = v0; v0 = v1; v1 = t; dv0 = dv1; dv1 = dt;         // v0,v1 = v1,v0 - q*v1
        //printf("r0 = "); polprint64s(r0,dr0);
        //printf("r1 = "); polprint64s(r1,dr1);
        //printf("v0 = "); polprint64s(v0,dv0);
        //printf("v1 = "); polprint64s(v1,dv1);
    }
    if( dv1>=0 ) {
        polcopy64s(v1,dv1,L);
        monic64s(L,dv1,p);
    }
    return dv1;
}

/*
Function: cppBM
What it does:
  C-compatible wrapper for BerlekampMassey64s used by Maple.
Inputs:
  - aIn/aLen: input sequence.
  - p: modulus.
  - lOut/outLen: output buffer.
  - degOut: output degree pointer.
Outputs:
  - Returns 0 on a valid call; negative codes report bad arguments or insufficient output space. degOut=-1 denotes no connection polynomial.
Example:
  int rc = cppBM(N,seq,101,N/2+1,L,&d);
*/
extern "C" int cppBM(int aLen,
    const LONG *aIn,
    const LONG p,
    int outLen,
    LONG *lOut,
    int *degOut)
{

// INITIAL CHECKS:

if (!aIn || !lOut || !degOut) {
    return -1;
}
if (aLen <= 0 || outLen <= 0) {
    return -2;
}

const int N = aLen;
const int n = N / 2;

// Lambda(x) has degree at most n = floor(N/2), so we need n+1 coefficient slots.
if (outLen < n + 1) {
    return -3;
}

*degOut = -1;
for (int i = 0; i < outLen; ++i) {
    lOut[i] = 0;
}

vector<LONG> a(aIn, aIn + N);
vector<LONG> Lbuf(N + 2, 0);
vector<LONG> W(3 * N + 8, 0);

int d = BerlekampMassey64s(a.data(), N, Lbuf.data(), W.data(), p);

if (d < 0) {
    // No connection polynomial (all-zero input or N==0 after stripping).
    // Not an error: report via degOut = -1 and return success.
    *degOut = -1;
    return 0;
}

if (d >= outLen) {
    *degOut = -1;
    return -4;
}

std::copy_n(Lbuf.data(), d + 1, lOut);
*degOut = d;
return 0;
}

/*
Function: newtonInterpMulRec
What it does:
  Interpolates the unique polynomial through n finite-field points using Newton divided differences, then converts the result to the monomial basis using mul64b.
Inputs:
  - x: n interpolation nodes.
  - y: n values; overwritten first by Newton coefficients and finally by monomial coefficients.
  - n: number of points.
  - p: modulus.
Outputs:
  - Returns the degree of the interpolating polynomial, or -1 when n<1 or a repeated node causes a zero denominator.
Example:
  int d = newtonInterpMulRec(x,y,n,101);  // y now stores coefficients
*/
int newtonInterpMulRec(LONG* x,
    LONG* y,
    const int n,
    const LONG p){
    if(n<1){
        return -1;
    }
    LONG *X=x;
    LONG *Y=y;
    int d;
    int i;
    int j;
    LONG prod;
    LONG s;
    for(j=1;j<n;j++){
        const LONG xj=X[j];
        s=Y[0];
        prod=sub64b(xj,X[0],p);
        for(i=1;i<j;i++){
            s=add64b(s,mul64b(prod,Y[i],p),p);
	        prod=mul64b(prod,sub64b(X[j],X[i],p),p);
        }
        if(prod==0){
            return -1;
        }
        Y[j]=mul64b(sub64b(Y[j],s,p),modinv64b(prod,p),p);	
    }
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            Y[j]=sub64b(Y[j],mul64b(X[d-i],Y[j+1],p),p);        
	    }
    }
    return d;
}

/*
Function: newtonInterpMulNormal
What it does:
  Newton interpolation kernel using the assembly modular multiplier on hot multiplication steps, then converting to monomial coefficients in place.
Inputs:
  - x: interpolation nodes.
  - y: values and output coefficient array.
  - n: number of points.
  - p: modulus.
Outputs:
  - Returns polynomial degree, or -1 on empty/repeated-node input; y is overwritten with monomial coefficients.
Example:
  int d = newtonInterpMulNormal(x,y,n,101);
*/
int newtonInterpMulNormal(LONG* x,
    LONG* y,
    const int n,
    const LONG p){
    if(n<1){
        return -1;
    }
    LONG *X=x;
    LONG *Y=y;
    int d;
    int i;
    int j;
    LONG prod;
    LONG s;
    for(j=1;j<n;j++){
        const LONG xj=X[j];
        s=Y[0];
        prod=sub64b(xj,X[0],p);
        for(i=1;i<j;i++){
            s=add64b(s,mul64bASM(prod,Y[i],p),p);
	        prod=mul64bASM(prod,sub64b(X[j],X[i],p),p);
        }
        if(prod==0){
            return -1;
        }
        Y[j]=mul64bASM(sub64b(Y[j],s,p),modinv64b(prod,p),p);	
    }
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            Y[j]=sub64b(Y[j],mul64bASM(X[d-i],Y[j+1],p),p);        
	    }
    }
    return d;
}

/*
Function: mkM
What it does:
  Builds the monic node polynomial M(x)=product_i (x-xs[i]) over GF(p) by repeated in-place multiplication.
Inputs:
  - m: destination/workspace vector; initialize m[0]=1 and allocate enough slots.
  - xs: interpolation nodes.
  - p: modulus.
Outputs:
  - Returns degree(M), with coefficients written into m.
Example:
  vector<LONG> M(xs.size()+1,0); M[0]=1; int d=mkM(M,xs,101);
*/
int mkM(vector<LONG>&m,const vector<LONG> &xs,const LONG p){    
    int degM=0;
    std::vector<LONG>linF(2,0);
    linF[1]=1;
    for(int i=0;i<xs.size();i++){
        linF[0]=(xs[i]==0?0:p-xs[i]);
        degM=pMULIP64(m.data(),linF.data(),degM,1,p);
    };
    return degM;
}


/*
NEWTON INTERPOLATION WRAPPER FOR MAPLE.
*/

/*
Function: cppInterp
What it does:
  C-compatible Maple wrapper for Newton interpolation.
Inputs:
  - xIn/xLen: interpolation nodes.
  - yIn/yLen: values.
  - p: modulus.
  - yOut/outLen: coefficient output buffer.
  - degOut: returned degree.
Outputs:
  - Returns 0 on success; negative codes report invalid dimensions, repeated-node failure, or insufficient space. yOut stores ascending monomial coefficients.
Example:
  int rc = cppInterp(n,x,n,y,101,n,coeff,&d);
*/
extern "C" int cppInterp(int xLen,
    const LONG *xIn,
    int yLen,
    const LONG *yIn,
    const LONG p,
    int outLen,
    LONG *yOut,
    int *degOut)
{


// INITIAL CHECKS: 

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

// INITIALIZING OUTPUT ARRAY: 

*degOut = -1;
for (int i = 0; i < outLen; ++i) {
    yOut[i] = 0;
}

// local working copies since kernel overwrites y
vector<LONG> x(xIn, xIn + n);
vector<LONG> y(yIn, yIn + n);

int d = newtonInterpMulNormal(x.data(), y.data(), n, p);

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

/*
Function: ftr_buildM
What it does:
  Internal helper that constructs M(x)=product_i (x-alpha[i]) for fault-tolerant reconstruction.
Inputs:
  - alpha: array of n nodes.
  - n: number of nodes.
  - p: modulus.
Outputs:
  - Returns a vector of n+1 coefficients for the monic degree-n product polynomial.
Example:
  auto M = ftr_buildM(alpha,n,101);
*/
static std::vector<LONG> ftr_buildM(const LONG *alpha,int n,LONG p){
    std::vector<LONG> M(n+1,0); M[0]=1; int degM=0; LONG lin[2]; lin[1]=1;
    for(int i=0;i<n;i++){ lin[0]=(alpha[i]==0?0:p-alpha[i]); degM=pMULIP64(M.data(),lin,degM,1,p); }
    return M;                                  // degree n
}

/*
Function: cppInterpDFTRFR
What it does:
  Combines Newton interpolation and deterministic fault-tolerant rational reconstruction in one C entry point.
Inputs:
  - alpha,Yin: nPts sample nodes and values.
  - N,D,E: reconstruction bounds.
  - p: modulus.
  - nOut,dOut: numerator and denominator output buffers.
Outputs:
  - Returns 0 on successful reconstruction, 1 when reconstruction fails, and a negative code for invalid input/interpolation failure.
Example:
  int rc = cppInterpDFTRFR(n,alpha,Y,N,D,E,101,N+1,num,D+1,den);
*/
extern "C" int cppInterpDFTRFR(int nPts,const LONG *alpha,const LONG *Yin,
                               int N,int D,int E,const LONG p,
                               int nOutLen,LONG *nOut,
                               int dOutLen,LONG *dOut){
    if(!alpha||!Yin||!nOut||!dOut)               return -1;
    if(nPts<1||nOutLen<N+1||dOutLen<D+1)         return -1;
    for(int i=0;i<nOutLen;i++) nOut[i]=0;
    for(int i=0;i<dOutLen;i++) dOut[i]=0;

    std::vector<LONG> xw(alpha,alpha+nPts), yw(Yin,Yin+nPts);
    int degU=newtonInterpMulNormal(xw.data(),yw.data(),nPts,p);   // U in yw[0..degU]
    if(degU<-1) return -6;                       // interpolation failure (repeated node)

    std::vector<LONG> M=ftr_buildM(alpha,nPts,p); int degM=nPts;
    std::vector<LONG> U;
    if(degU>=0) U.assign(yw.begin(),yw.begin()+degU+1); else U.assign(1,0);

    int degNum=-1,degDen=-1;
    int st=dftrfr(M,U,degM,degU,N,D,E,p,nOut,degNum,dOut,degDen);
    return (st==0)?1:0;
}

/*
Function: cppAffineLine
What it does:
  Generates T parameter points on the affine line used by MRFI, using alpha as the first coordinate and beta/sigma for the remaining coordinates.
Inputs:
  - T,numVar: number of points and parameters.
  - alpha,beta,sigma: line-defining arrays with supplied lengths.
  - p: modulus.
  - out/outLen: flat point-major destination array.
Outputs:
  - Returns 0 on success or -1 for invalid arguments; writes T*numVar residues to out.
Example:
  int rc = cppAffineLine(T,nv,T,alpha,nv-1,beta,nv,sigma,101,T*nv,out);
*/
extern "C" int cppAffineLine(int T,int numVar,
                             int alphaLen,const LONG *alpha,
                             int betaLen,const LONG *beta,
                             int sigmaLen,const LONG *sigma,
                             const LONG p,
                             int outLen,LONG *out){
    if(!alpha||!beta||!sigma||!out)                       return -1;
    if(T<1||numVar<1)                                     return -1;
    if(alphaLen<T)                                        return -1;
    if(numVar>1 && (betaLen<numVar-1||sigmaLen<numVar))   return -1;
    if((long long)outLen < (long long)T*numVar)           return -1;
    LONG s1=sigma[0]%p; if(s1<0) s1+=p;
    for(int np=0;np<T;np++){
        LONG a=alpha[np]%p; if(a<0) a+=p;
        LONG d=sub64b(a,s1,p);                 /* alpha[np] - sigma[1] mod p */
        LONG *row=out+(size_t)np*numVar;
        row[0]=a;
        for(int j=1;j<numVar;j++){
            LONG b=beta[j-1]%p;  if(b<0) b+=p;
            LONG c=sigma[j]%p;   if(c<0) c+=p;
            row[j]=add64b(mul64b(b,d,p),c,p);
        }
    }
    return 0;
}

/*
Function: cppFTREval
What it does:
  Interpolates sample values, performs deterministic fault-tolerant rational reconstruction, and evaluates the recovered numerator and denominator at sigma.
Inputs:
  - alpha,Yin: sample nodes/values and their capacities.
  - nPts: number of active samples.
  - sigma: evaluation point.
  - N,D,E: reconstruction bounds.
  - p: modulus.
  - out: at least two slots.
Outputs:
  - Returns 0 on success, 1 if reconstruction fails, or a negative code on invalid/repeated-node input; out[0]=num(sigma), out[1]=den(sigma).
Example:
  LONG out[2]; int rc=cppFTREval(n,n,alpha,n,Y,sigma,N,D,E,101,2,out);
*/
extern "C" int cppFTREval(int nPts,
                          int alphaLen,const LONG *alpha,
                          int yLen,const LONG *Yin,
                          LONG sigma,int N,int D,int E,const LONG p,
                          int outLen,LONG *out){
    if(!alpha||!Yin||!out)        return -1;
    if(nPts<1||outLen<2)          return -1;
    if(nPts>alphaLen||nPts>yLen)  return -1;   /* buffers may be larger than nPts */
    if(N<0||D<0||E<0)             return -1;
    out[0]=0; out[1]=0;

    std::vector<LONG> xw(alpha,alpha+nPts), yw(Yin,Yin+nPts);
    int degU=newtonInterpMulNormal(xw.data(),yw.data(),nPts,p);
    if(degU<-1) return -6;                       /* repeated node */

    std::vector<LONG> M=ftr_buildM(alpha,nPts,p); int degM=nPts;
    std::vector<LONG> U;
    if(degU>=0) U.assign(yw.begin(),yw.begin()+degU+1); else U.assign(1,0);

    std::vector<LONG> num(N+1,0),den(D+1,0);
    int degNum=-1,degDen=-1;
    if(dftrfr(M,U,degM,degU,N,D,E,p,num.data(),degNum,den.data(),degDen)==0)
        return 1;                                /* reconstruction failed */

    out[0]=pEVAL64(num.data(),degNum,sigma,p);
    out[1]=pEVAL64(den.data(),degDen,sigma,p);
    return 0;
}

/*
Function: ratRECON_C
What it does:
  C-compatible rational-reconstruction wrapper around ratReconNormal, used by the Maple mRATRECON binding.
Inputs:
  - M/U arrays with lengths/degrees.
  - N,D: requested degree bounds.
  - p: modulus.
  - nOut,dOut: coefficient buffers.
  - degNOUT,degDOUT: returned degrees.
Outputs:
  - Returns 0 on success or a negative error code; writes normalized numerator and denominator coefficients and their degrees.
Example:
  int rc = ratRECON_C(mLen,dM,M,uLen,dU,U,N,D,101,nLen,num,&dn,dLen,den,&dd);
*/
extern "C" int ratRECON_C(int mLen,
    int degM,
    const LONG *M,
    int uLen,
    int degU,
    const LONG *U,
    const int N,
    const int D,
    const LONG p,
    int nOutLen,
    LONG *nOut,
    int *degNOUT,
    int dOutLen,
    LONG *dOut,
    int *degDOUT)
{

if (!M || !U || !nOut || !dOut || !degNOUT || !degDOUT) return -1;
if (degM < 0 || degU < 0) return -1;
if (mLen <= 0 || uLen <= 0 || nOutLen <= 0 || dOutLen <= 0) return -1;
if (degM >= mLen || degU >= uLen) return -1;

*degNOUT = -1;
*degDOUT = -1;

for (int i = 0; i < nOutLen; ++i) nOut[i] = 0;
for (int i = 0; i < dOutLen; ++i) dOut[i] = 0;

int wsSize = std::max(degM, degU) + 1;

std::vector<LONG> m(M, M + (degM + 1));
std::vector<LONG> u(U, U + (degU + 1));
std::vector<LONG> rTmp(wsSize, 0);
std::vector<LONG> tTmp(wsSize, 0);

int degROut = -1;
int degTOut = -1;
RatReconFastWS W(wsSize);

int rc = ratReconNormal(m,
            u,
            degM,
            degU,
            N,
            D,
            p,
            W,
            rTmp.data(),
            degROut,
            tTmp.data(),
            degTOut);

if (rc != 0) {
*degNOUT = -1;
*degDOUT = -1;
return rc;
}

if (degROut < 0 || degTOut < 0) return -2;
if (degROut >= wsSize || degTOut >= wsSize) return -3;
if (degROut + 1 > nOutLen) return -4;
if (degTOut + 1 > dOutLen) return -5;
std::copy_n(rTmp.data(), degROut + 1, nOut);
std::copy_n(tTmp.data(), degTOut + 1, dOut);
*degNOUT = degROut;
*degDOUT = degTOut;

return 0;
}

/* 
Mike's linear Algebra routines for fast computation.
*/

/*
Function: rref
What it does:
  Reduces an n-by-m row-major matrix to reduced row echelon form over GF(p).
Inputs:
  - B: flat matrix array, overwritten in place; entries are assumed in 0..p-1.
  - n,m: row and column counts.
  - p: modulus.
Outputs:
  - Returns rank(B); B contains its RREF on return.
Example:
  LONG rank = rref(B,n,m,101);
*/
LONG rref( LONG *B, int n, int m, LONG p ) {
// Put B in reduced row Echelon form and return rank(B)
// The code assumes 0 <= B[i,j] < p
LONG t,det;
int c,r,i,j;
//recint P;
   //printf("rref: n=%d m=%d\n",n,m);
   //P = recip1(p);
   det = 1;
   for( c=0,r=0; c<m && r<n; c++ ) {
      // Search for a pivot element
      for( i=r; i<n && B[m*i+c]==0; i++ );
      if( i==n ) { det = 0; continue; }
      if( i!=r ) { // interchange row i with row r
          det = neg64s(det,p);
          for( j=c; j<m; j++ ) { t = B[i*m+j]; B[i*m+j] = B[r*m+j]; B[r*m+j] = t; }
      }
      det = mul64b(det,B[r*m+c],p);
      t = modinv64b(B[r*m+c],p);
      for( j=c+1; j<m; j++ ) B[r*m+j] = mul64b(t,B[r*m+j],p);
      B[r*m+c] = 1;
      for( i=0; i<n; i++ ) {
         if( i==r || B[i*m+c]==0 ) continue;
         for( j=c+1; j<m; j++ )
            B[i*m+j] = sub64b(B[i*m+j],mul64b(B[i*m+c],B[r*m+j],p),p);
         B[i*m+c] = 0;
      }
      r++;  // go to next row
   }
   //printf("det(B)=%lld\n",det);
   return(r); // r = rank(B)
}

/*
RREF WRAPPER FOR MAPLE.

B is the n x m matrix in C (row major) order.  It is OVERWRITTEN with its
reduced row echelon form.  The return value is rank(B), or a negative code
on bad arguments.  Entries are reduced into [0,p) first, since rref assumes
0 <= B[i,j] < p.
*/

/*
Function: cppRREF
What it does:
  C-compatible wrapper for rref that validates dimensions and first normalizes all matrix entries into 0..p-1.
Inputs:
  - n,m: matrix dimensions.
  - B: row-major matrix overwritten in place.
  - p: modulus.
Outputs:
  - Returns the rank on success or a negative argument error code.
Example:
  int rank = cppRREF(n,m,B,101);
*/
extern "C" int cppRREF(int n,
    int m,
    LONG *B,
    const LONG p)
{

// INITIAL CHECKS:

if (!B) {
    return -1;
}
if (n <= 0 || m <= 0) {
    return -2;
}
if (p < 2) {
    return -3;
}

// NORMALIZE INPUT INTO [0,p):

for (size_t i = 0; i < (size_t)n*m; ++i) {
    LONG v = B[i] % p;
    B[i] = v + ((v>>63)&p);
}

return (int)rref(B, n, m, p);
}

/*
PARAMETRIC MATRIX EVALUATION + SOLVE FOR MAPLE.

The black box evaluates one fixed augmented matrix L(y1,...,ynv) at a point
and solves the resulting system mod p, millions of times.  Doing the
evaluation in Maple costs one interpreted eval + one modp per entry, which
measured at about 1.5 microseconds per entry, i.e. more than the whole
linear solve.  Here L is encoded ONCE in Maple as a sparse list of monomials
and the per call work is a gather plus nv multiplications per term.

ENCODING (all arrays flat, row major, entry (i,j) has index i*nc+j):

  entStart[0..nr*nc]   CSR style: the terms of entry e are
                       t = entStart[e] .. entStart[e+1]-1
  coef[t]              coefficient of term t, already reduced into [0,p)
  expo[t*nv+k]         exponent of parameter k in term t
  pnt[k]               value of parameter k
  dmax                 max exponent appearing in expo

Entry (i,j) evaluates to  sum_t coef[t] * prod_k pnt[k]^expo[t*nv+k] mod p.
The powers are tabulated once per call, so a term costs nv multiplies and
no exponentiation.

Return value:  0  unique solution written to xOut
               1  singular or inconsistent  (Maple side returns FAIL)
              <0  bad arguments
info[0] = rank of the evaluated augmented matrix.
*/

/*
Shared kernel: evaluate the encoded matrix at ONE point and solve.
B and pw are caller supplied workspaces so the batch entry point below can
allocate them once for the whole block.
  returns  0 unique solution written to x
           1 singular or inconsistent
  rank is written to *rankOut.
*/

/*
Function: evalSolveOne
What it does:
  Internal black-box kernel that evaluates a sparsely encoded parametric augmented matrix at one parameter point and solves it by RREF.
Inputs:
  - nr,nc,nv: row, augmented-column, and parameter counts.
  - entStart,expo,coef: sparse monomial encoding.
  - pnt: parameter values.
  - p,dmax: modulus and maximum exponent.
  - x,rankOut: solution and rank outputs.
  - pw,B: caller-supplied power-table and matrix workspaces.
Outputs:
  - Returns 0 for a unique solution and 1 for a singular/inconsistent system; writes x and rankOut.
Example:
  int st = evalSolveOne(nr,nr+1,nv,ent,expo,coef,pnt,p,dmax,x,&rank,pw,B);
*/
static int evalSolveOne(int nr,int nc,int nv,
    const int *entStart,const int *expo,const LONG *coef,
    const LONG *pnt,const LONG p,int dmax,
    LONG *x,int *rankOut,LONG *pw,LONG *B)
{
const int nEnt = nr*nc;
const int W = dmax+1;

// POWER TABLE: pw[k*W+d] = pnt[k]^d mod p

for (int k = 0; k < nv; ++k) {
    LONG a = pnt[k];
    if (a < 0 || a >= p) { a %= p; a += (a>>63)&p; }
    LONG *row = pw + (size_t)k*W;
    row[0] = 1;
    for (int d = 1; d < W; ++d) row[d] = mul64b(row[d-1],a,p);
}

// EVALUATE THE MATRIX.  The coefficients arrive already reduced into [0,p),
// so the test below costs a predicted branch instead of a divq per term.

for (int e = 0; e < nEnt; ++e) {
    LONG acc = 0;
    const int t0 = entStart[e];
    const int t1 = entStart[e+1];
    for (int t = t0; t < t1; ++t) {
        LONG c = coef[t];
        if (c < 0 || c >= p) { c %= p; c += (c>>63)&p; }
        const int *ex = expo + (size_t)t*nv;
        for (int k = 0; k < nv && c != 0; ++k) {
            const int d = ex[k];
            if (d) c = mul64b(c,pw[(size_t)k*W+d],p);
        }
        acc = add64b(acc,c,p);
    }
    B[e] = acc;
}

// SOLVE

int r = (int)rref(B, nr, nc, p);
if (rankOut) *rankOut = r;
if (r != nr) {
    return 1;                                  // rank deficient
}
for (int i = 0; i < nr; ++i) {
    if (B[(size_t)i*nc+i] != 1) return 1;      // pivot in the b column
}
for (int i = 0; i < nr; ++i) {
    x[i] = B[(size_t)i*nc+nr];
}
return 0;
}

/*
Function: cppEvalSolve
What it does:
  C-compatible single-point parametric black-box wrapper around evalSolveOne with reusable static workspaces.
Inputs:
  - nr,nc,nv and sparse encoding arrays.
  - pnt: one parameter point.
  - p,dmax: modulus and maximum exponent.
  - xOut/outLen: solution buffer.
  - info/infoLen: metadata buffer.
Outputs:
  - Returns 0 for a unique solution, 1 for singular/inconsistent input, or a negative validation code; info[0] receives the rank.
Example:
  int rc = cppEvalSolve(nr,nr+1,nv,ent,expo,coef,pnt,p,dmax,nr,x,1,info);
*/
extern "C" int cppEvalSolve(int nr,
    int nc,
    int nv,
    const int *entStart,
    const int *expo,
    const LONG *coef,
    const LONG *pnt,
    const LONG p,
    int dmax,
    int outLen,
    LONG *xOut,
    int infoLen,
    LONG *info)
{

// INITIAL CHECKS:

if (!entStart || !pnt || !xOut || !info) {
    return -1;
}
if (nr <= 0 || nc != nr+1 || nv <= 0 || dmax < 0) {
    return -2;
}
if (outLen < nr || infoLen < 1) {
    return -3;
}
if (p < 2) {
    return -4;
}

const int nEnt = nr*nc;
const int nT   = entStart[nEnt];
if (nT < 0) {
    return -5;
}
if (nT > 0 && (!expo || !coef)) {
    return -1;                             // coef/expo may be NULL only if L == 0
}

// INITIALIZING OUTPUT ARRAY:

for (int i = 0; i < outLen; ++i) {
    xOut[i] = 0;
}
info[0] = 0;

// Workspaces are static so that the black box does not malloc on every
// call.  Maple runs single threaded here (kernelopts(numcpus=1)).

static std::vector<LONG> pw;
static std::vector<LONG> B;

const int W = dmax+1;
if ((int)pw.size() < nv*W) pw.resize((size_t)nv*W);
if ((int)B.size() < nEnt)  B.resize((size_t)nEnt);

int rank = 0;
int st = evalSolveOne(nr,nc,nv,entStart,expo,coef,pnt,p,dmax,
                      xOut,&rank,pw.data(),B.data());
info[0] = rank;
return st;
}

/*
BATCHED, TRANSPOSED BLACK BOX.

The MRFI loop asks for the whole block of npts points on an affine line and
then, for each equation i, needs the sequence of x_i over those points.  Done
one point at a time in Maple that is npts external calls, npts solution lists,
and an interpreted transpose of an npts x nr matrix.  Here the block is solved
in ONE call and the answers are written already transposed:

    xOut[(i-1)*ldx + (s-1)]  =  x_i evaluated at point s

so row i is contiguous and Maple can hand it straight to cppFTREval through
ArrayTools:-Alias with no copy at all.

pts holds the points as cppAffineLine writes them, point major:
    pts[(s-1)*nv + k]  =  parameter k+1 of point s

Return value:  0  every point solved
               1  at least one point was singular or inconsistent
              <0  bad arguments
info[0] = number of points solved, info[1] = 1-based index of the first
singular point, or 0 when there was none.
*/

/*
Function: cppEvalSolveBlock
What it does:
  Evaluates and solves the same encoded augmented system at a block of points, storing solutions already transposed by equation.
Inputs:
  - nr,nc,nv and sparse encoding arrays.
  - pts/npts: flat point-major parameter block.
  - p,dmax: modulus and maximum exponent.
  - ldx,xOut: row stride and transposed output buffer.
  - info: solved-count and first-failure metadata.
Outputs:
  - Returns 0 if every point solves, 1 if any point is singular/inconsistent, or a negative validation code. info[0]=solved count, info[1]=first bad 1-based point.
Example:
  int rc = cppEvalSolveBlock(nr,nr+1,nv,ent,expo,coef,pts,npts,p,dmax,npts,X,2,info);
*/
extern "C" int cppEvalSolveBlock(int nr,
    int nc,
    int nv,
    const int *entStart,
    const int *expo,
    const LONG *coef,
    const LONG *pts,
    int npts,
    const LONG p,
    int dmax,
    int ldx,
    LONG *xOut,
    int infoLen,
    LONG *info)
{

// INITIAL CHECKS:

if (!entStart || !pts || !xOut || !info) {
    return -1;
}
if (nr <= 0 || nc != nr+1 || nv <= 0 || dmax < 0 || npts <= 0) {
    return -2;
}
if (ldx < npts || infoLen < 2) {
    return -3;
}
if (p < 2) {
    return -4;
}

const int nEnt = nr*nc;
const int nT   = entStart[nEnt];
if (nT < 0) {
    return -5;
}
if (nT > 0 && (!expo || !coef)) {
    return -1;
}

info[0] = 0;
info[1] = 0;

static std::vector<LONG> pw;
static std::vector<LONG> B;
static std::vector<LONG> x;

const int W = dmax+1;
if ((int)pw.size() < nv*W) pw.resize((size_t)nv*W);
if ((int)B.size() < nEnt)  B.resize((size_t)nEnt);
if ((int)x.size() < nr)    x.resize((size_t)nr);

int solved = 0;
int firstBad = 0;

for (int s = 0; s < npts; ++s) {
    int rank = 0;
    int st = evalSolveOne(nr,nc,nv,entStart,expo,coef,
                          pts + (size_t)s*nv,p,dmax,
                          x.data(),&rank,pw.data(),B.data());
    if (st != 0) {
        if (firstBad == 0) firstBad = s+1;
        for (int i = 0; i < nr; ++i) xOut[(size_t)i*ldx + s] = 0;
        continue;
    }
    for (int i = 0; i < nr; ++i) xOut[(size_t)i*ldx + s] = x[i];
    solved++;
}

info[0] = solved;
info[1] = firstBad;

return firstBad ? 1 : 0;
}


/*
COMPATIBILITY LAYER FOR THE ROOT FINDING CODE.

polgcdext64s, polpowmod64s, polsplit64s and polroots64s are written against
Mike's ...64s naming.  Everything they need already exists in this file under
the ...64b / ...IP64 names, except for four routines that were not here at all.
Nothing below duplicates existing arithmetic: the forwarders compile away, and
only polgcd64s, polscamul64s, polsqr64s and polprint64s add code.

MISSING BEFORE, DEFINED HERE:
  max32s        integer max
  add64s sub64s mul64s modinv64s     names for add64b sub64b mul64b modinv64b
  poldiv64s     name for polDIVIP64  (same contract: a := [remainder,quotient],
                returns deg of the remainder, quotient starts at a+degB)

polgcdext64s used polsubmulP with a recint reciprocal.  It now calls
polSUBMUL64 directly, which does the same A -= (a*x+b)*B with ZMUL/ZFMA/ZMOD
accumulators and no reciprocal, so there is no recint anywhere in the root
finding path.
  polscamul64s  in place scalar multiply on a raw pointer
  polsqr64s     square, via polmul64s with both operands aliased
  polgcd64s     monic gcd, RESULT LEFT IN THE FIRST ARGUMENT
  polprint64s   only referenced from commented out debugging lines, but the
                calls are there, so it is defined rather than deleted
  ULONG         polroots64s declares  extern ULONG seed,mult;  and the globals
                at the top of this file are ULNG, so the two names must denote
                the same type or that declaration is a type mismatch
*/

using ULONG = ULNG;

/*
Function: max32s
What it does:
  Returns the larger of two int values; retained as a compatibility helper for the root-finding code.
Inputs:
  - a,b: integer values.
Outputs:
  - Returns max(a,b).
Example:
  int m = max32s(3,7);  // m = 7
*/
inline int max32s(int a,int b){ return a>b ? a : b; }

/*
Function: add64s
What it does:
  Compatibility alias that forwards modular addition to add64b.
Inputs:
  - a,b: residues.
  - p: modulus.
Outputs:
  - Returns (a+b) mod p.
Example:
  LONG c = add64s(70,50,101);
*/
inline LONG add64s(LONG a,LONG b,LONG p){ return add64b(a,b,p); }
/*
Function: sub64s
What it does:
  Compatibility alias that forwards modular subtraction to sub64b.
Inputs:
  - a,b: residues.
  - p: modulus.
Outputs:
  - Returns (a-b) mod p.
Example:
  LONG c = sub64s(3,5,101);
*/
inline LONG sub64s(LONG a,LONG b,LONG p){ return sub64b(a,b,p); }
/*
Function: mul64s
What it does:
  Compatibility alias that forwards modular multiplication to mul64b.
Inputs:
  - a,b: residues.
  - p: modulus.
Outputs:
  - Returns (a*b) mod p.
Example:
  LONG c = mul64s(25,9,101);
*/
inline LONG mul64s(LONG a,LONG b,LONG p){ return mul64b(a,b,p); }
/*
Function: modinv64s
What it does:
  Compatibility alias that forwards modular inversion to modinv64b.
Inputs:
  - c: value to invert.
  - p: modulus.
Outputs:
  - Returns c^(-1) mod p or 0 when no inverse exists.
Example:
  LONG inv = modinv64s(3,101);
*/
inline LONG modinv64s(LONG c,LONG p){ return modinv64b(c,p); }

/*
Function: poldiv64s
What it does:
  Compatibility alias for polDIVIP64 used by the root-finding routines.
Inputs:
  - a: dividend/workspace.
  - b: divisor.
  - da,db: degrees.
  - p: modulus.
Outputs:
  - Returns the remainder degree while leaving remainder/quotient in a's standard in-place division layout.
Example:
  int dr = poldiv64s(A,B,da,db,101);
*/
inline int poldiv64s(LONG *a,const LONG *b,int da,int db,LONG p){
    return polDIVIP64(a,b,da,db,p);
}

// A := c*A
/*
Function: polscamul64s
What it does:
  Scales a raw polynomial A by the field element c in place.
Inputs:
  - c: scalar.
  - A: coefficient array.
  - d: degree of A.
  - p: modulus.
Outputs:
  - Returns nothing; overwrites A with c*A mod p.
Example:
  polscamul64s(5,A,d,101);
*/
void polscamul64s(LONG c,LONG *A,int d,LONG p){
    int i;
    if( c==1 ) return;
    for( i=0; i<=d; i++ ) A[i]=mul64b(A[i],c,p);
    return;
}

// C := A^2.  polmul64s reads A twice and writes a separate C, so aliasing the
// two inputs is safe.  C must hold 2*da+1 coefficients.
/*
Function: polsqr64s
What it does:
  Squares a polynomial by calling polmul64s with the same source as both factors.
Inputs:
  - A: input coefficient array.
  - C: separate destination with at least 2*da+1 slots.
  - da: degree of A.
  - p: modulus.
Outputs:
  - Returns degree(C), where C=A^2 mod p.
Example:
  int dc = polsqr64s(A,C,da,101);
*/
int polsqr64s(LONG *A,LONG *C,int da,LONG p){
    return polmul64s(A,A,C,da,da,p);
}

/*
Function: polprint64s
What it does:
  Prints a raw coefficient-array polynomial for debugging.
Inputs:
  - A: coefficient array.
  - d: degree; negative prints 0.
Outputs:
  - Returns nothing; writes the polynomial to standard output.
Example:
  polprint64s(A,d);
*/
void polprint64s(LONG *A,int d){
    int i;
    if( d<0 ) { printf("0\n"); return; }
    for( i=d; i>=0; i-- ) {
        if( A[i]==0 ) continue;
        printf("%lld",(long long)A[i]);
        if( i>0 ) printf("*x^%d",i);
        if( i>0 ) printf(" + ");
    }
    printf("\n");
    return;
}

/*
Monic gcd(A,B) in Zp[x].  BOTH inputs are destroyed and the answer is left in
A, which is what polsplit64s and polroots64s assume:

    dg = polgcd64s( W, W+d, da, d, p );   // g = gcd(W,f) lands in W
    dg = polgcd64s( f, W, d, da, p );     // g = gcd(f,W-1) lands in f

Returns deg(gcd), so 0 means the gcd is the constant 1.
*/
/*
Function: polgcd64s
What it does:
  Computes the monic gcd of two raw polynomials over GF(p) for the root-finding path; both inputs may be destroyed.
Inputs:
  - A,B: mutable polynomial arrays.
  - da,db: degrees.
  - p: modulus.
Outputs:
  - Returns degree(gcd) and guarantees the monic gcd is stored in A.
Example:
  int dg = polgcd64s(A,B,da,db,101);
*/
int polgcd64s(LONG *A,LONG *B,int da,int db,LONG p){
    LONG *r1,*r2,*t;
    int d1,d2,dr,i;
    LONG inv;
    while( da>=0 && A[da]==0 ) da--;
    while( db>=0 && B[db]==0 ) db--;
    if( da<0 && db<0 ) return -1;                 // gcd(0,0)
    if( da<0 ) { polcopy64s(B,db,A); da=db; db=-1; }
    r1=A; d1=da; r2=B; d2=db;
    if( d1<d2 ) { t=r1; r1=r2; r2=t; i=d1; d1=d2; d2=i; }
    while( d2>=0 ) {
        dr = polDIVIP64(r1,r2,d1,d2,p);           // remainder in r1[0..dr]
        t=r1; r1=r2; r2=t;
        d1=d2; d2=dr;
    }
    if( d1<0 ) return -1;
    if( r1!=A ) polcopy64s(r1,d1,A);              // answer must end up in A
    if( A[d1]!=1 ) {
        inv=modinv64b(A[d1],p);
        for( i=0; i<d1; i++ ) A[i]=mul64b(A[i],inv,p);
        A[d1]=1;
    }
    return d1;
}

/*
Function: polgcdext64s
What it does:
  Computes an extended polynomial gcd over GF(p), producing G,S,T such that S*A + T*B = G with G monic.
Inputs:
  - A,B and da,db: nonzero inputs used destructively as remainder storage.
  - G: gcd output.
  - S,T: optional Bezout-coefficient outputs; pass null to omit either.
  - dG,dS,dT: returned degrees.
  - W: scratch workspace.
  - p: modulus.
Outputs:
  - Returns nothing; fills G and requested Bezout coefficients/degrees.
Example:
  polgcdext64s(A,B,da,db,G,S,T,&dG,&dS,&dT,W,101);
*/
void polgcdext64s( LONG *A, LONG *B, int da, int db,
    LONG *G, LONG *S, LONG *T, int *dG, int *dS, int *dT,
    //LONG *s1, *s2, *t1, *t2, int *ds1, int *ds2, int *dt1, int *dt2,
    LONG *W, LONG p )
{
// Solve S A + T B = G = monic gcd(A,B) for G,S,T in Zp[x]
// The arrays A and B are used for the remainder sequence so they are destroyed
// G,S,T must all be of size max(da+1,db+1)
// if( S==0 ) W is working storage of size max(da+1,db+1)
// if( T==0 ) W is working storage of size max(da+1,db+1)
// if S==0 or T==0 then S (and/or T) are not computed

int m,dr,ds,dt,dq,ds1,ds2,dt1,dt2; LONG a,b,u;
LONG *q,*r,*r1,*r2,*s,*s1,*s2,*t,*t1,*t2;

if( da<0 || db<0 ) { printf("inputs must be non-zero\n"); exit(1); }
m = max32s(da+1,db+1);
r1 = A; r2 = B;
if(S) { s1 = S; s2 = W;   s1[0]=1; ds1=0; ds2=-1; }
if(T) { t1 = T; 
if(S) t2 = W+m; else t2 = W;
t2[0]=1; dt2=0; dt1=-1;
}
while( 1 ) {
if( db>0 && da-db==1 ) { // normal case
u = modinv64s(r2[db],p);
a = mul64s(r1[da],u,p);
b = mul64s(a,r2[db-1],p);
b = mul64s(u,sub64s(r1[da-1],b,p),p);             // quotient = a x + b
//dr = polsubmul(r1,r2,a,b,da,db,p);                // r1 = r1 - (a x + b) r2
//if(S) ds = polsubmul(s1,s2,a,b,ds1,ds2,p);        // s1 = s1 - (a x + b) s2
//if(T) dt = polsubmul(t1,t2,a,b,dt1,dt2,p);        // t1 = t1 - (a x + b) t2
dr = polSUBMUL64(r1,r2,a,b,da,db,p);               // r1 = r1 - (a x + b) r2
if(S) ds = polSUBMUL64(s1,s2,a,b,ds1,ds2,p);       // s1 = s1 - (a x + b) s2
if(T) dt = polSUBMUL64(t1,t2,a,b,dt1,dt2,p);       // t1 = t1 - (a x + b) t2
}
else {
dr = poldiv64s(r1,r2,da,db,p);                 // r1 = [remainder,quotient]
q  = r1+db; dq = da-db;
if(S) ds = polmul64s(q,s2,G,dq,ds2,p);
if(S) ds = polsub64s(s1,G,s1,ds1,ds,p);        // s1 = s1 - q s2
if(T) dt = polmul64s(q,t2,G,dq,dt2,p);
if(T) dt = polsub64s(t1,G,t1,dt1,dt,p);        // t1 = t1 - q t2
}
if( dr<0 ) { /* D|C so gcd(A,B)=D */
polcopy64s(r2,db,G);
if(S) if( s2!=S ) polcopy64s(s2,ds2,S);
if(T) if( t2!=T ) polcopy64s(t2,dt2,T);
if( G[db]!=1 ) {
  u = modinv64s(G[db],p);
  polscamul64s(u,G,db,p);
  if(S) polscamul64s(u,S,ds2,p);
  if(T) polscamul64s(u,T,dt2,p);
}
dG[0] = db;
if(S) dS[0] = ds2;
if(T) dT[0] = dt2;
return;
}
r = r1; r1 = r2; r2 = r;  da = db;   db = dr;
if(S) { s = s1; s1 = s2; s2 = s; ds1 = ds2; ds2 = ds; }
if(T) { t = t1; t1 = t2; t2 = t; dt1 = dt2; dt2 = dt; }
}
}

// S = 1/A mod B
/*
Function: polmodinv64s
What it does:
  Computes the polynomial inverse S = 1/A mod M when A is invertible modulo M, using the extended gcd.
Inputs:
  - A,M: input polynomials.
  - da,dm: degrees with da<dm.
  - G,S: gcd and inverse/cofactor outputs.
  - W: scratch storage for copied inputs and EEA work.
  - p: modulus.
Outputs:
  - Returns degree(S) when gcd(A,M)=1; returns -deg(G) when the gcd has positive degree.
Example:
  int ds = polmodinv64s(A,M,da,dm,G,S,W,101);
*/
int polmodinv64s( LONG *A, LONG *M, int da, int dm,
    LONG *G, LONG *S, LONG *W, LONG p )
{   int dG, dS, dT;
// W must be able to hold three polynomials of degree dm
while( da>=0 && A[da]==0 ) da--;
if( da<0 ) { printf("A is zero\n"); exit(1); }
while( dm>=0 && M[dm]==0 ) dm--;
if( da>=dm ) { printf("deg(A) < deg(M) error\n"); exit(1); }
//printf("dA=%d dM=%d\n",da,dm);
// copy A and M into W
polcopy64s(A,da,W); A = W; W += da+1;
polcopy64s(M,dm,W); M = W; W += dm+1;
polgcdext64s( A, M, da, dm, G, S, 0, &dG, &dS, &dT, W, p );
//printf("dG=%d  dS=%d\n",dG,dS);
if( dG>0 ) return -dG; else return dS;
}

/* C(x) := A(x)^n mod B(x) mod p;  0<=deg(A)<deg(B) and R must be of size 2*db-1 */
/* If A(x) is not reduced mod B(x) then we first compute C(x) := A(x) mod B(x)   */
/*
Function: polpowmod64s
What it does:
  Computes C(x)=A(x)^n mod B(x) over GF(p) by repeated squaring and modular polynomial reduction.
Inputs:
  - A: base polynomial, possibly reduced in place if degree >= db.
  - n: nonnegative exponent.
  - B: modulus polynomial.
  - da,db: degrees.
  - C: output array.
  - R: scratch array large enough for intermediate products.
  - p: modulus.
Outputs:
  - Returns degree(C) and writes the reduced power into C.
Example:
  int dc = polpowmod64s(A,10,B,da,db,C,R,101);
*/
int polpowmod64s( LONG * A, LONG n, LONG * B, int da, int db, LONG *C, LONG *R, LONG p )
{
    int dc,k,b[63];

    if( n==0 ) { C[0] = 1; return 0; }
    if( da>=db ) da = poldiv64s(A,B,da,db,p);                   // reduce A mod B first
    for( k=0; n>0; k++ ) { b[k]=n&1; n=n/2; }
    polcopy64s(A,da,C);
    dc = da;
    k--;
    while( k>0 ) { k--;
       // Main step: compute C := C^2 mod B in Zp[x]
       //dc = polmul64s(C,C,R,dc,dc,p);                           //printf("deg(R) = %d; R = ",dc); polprint64s(R,dc);
       dc = polsqr64s(C,R,dc,p);                                //printf("deg(R) = %d; R = ",dc); polprint64s(R,dc);
       dc = poldiv64s(R,B,dc,db,p);
       polcopy64s(R,dc,C);                                      //printf("deg(C) = %d; C = ",dc); polprint64s(C,dc);
       if( b[k]==1 ) {                                          //printf(" b[%d]=%d \n", k, b[k] );
           dc = polmul64s(A,C,R,da,dc,p);                       //printf("deg(R) = %d; R = ",dc); polprint64s(R,dc);
           dc = poldiv64s(R,B,dc,db,p);
           polcopy64s(R,dc,C);                                  //printf("deg(C) = %d; C = ",dc); polprint64s(C,dc);
       }
    }
    return dc;
}

// Input f in Zp[x] of degree d > 0, a known product of d linear factors.
// Output roots of f in R.
// The input array f is destroyed.
// W is a scratch array of size at least 3*d
/*
Function: polsplit64s
What it does:
  Recursively splits a polynomial known to be a product of distinct linear factors over GF(p), using randomized gcd separation.
Inputs:
  - f: degree-d polynomial, destroyed during splitting.
  - d: positive degree.
  - R: output array for d roots.
  - W: scratch array of at least 3*d elements.
  - p: modulus.
Outputs:
  - Returns nothing; writes the d roots into R.
Example:
  polsplit64s(f,d,roots,W,101);
*/
void polsplit64s( LONG *f, int d, LONG *R, LONG *W, LONG p )
{
   int da,dg; LONG alpha, A[2];
   if( d==1 ) { alpha = p-f[0]; R[0] = alpha; return; }
   alpha = rand64s(p); A[1] = 1; A[0] = alpha;
   da = polpowmod64s( A, (p-1)/2, f, 1, d, W, W+d, p );
   if( da==0 ) return polsplit64s(f,d,R,W,p);      // alpha is unlucky, try again
   W[0] = add64s(W[0],1,p);                        // W = (x+alpha)^((p-1)/2) + 1 mod f
   polcopy64s( f, d, W+d );
   dg = polgcd64s( W, W+d, da, d, p );             // g = gcd( W, f ) in W
   if( dg==0 ) return polsplit64s(f,d,R,W,p);      // g = 1 ==> alpha is unlucky, try again
   poldiv64s(f,W,d,dg,p);                          // compute quotient q = f/g destroying f
   polcopy64s(W,dg-1,f);                           // f = [ g mod x^dg followed by q ]
   polsplit64s(f+dg,d-dg,R,W,p);
   f[dg] = 1;
   polsplit64s(f,dg,R+d-dg,W,p);
   return;
}

/*
Function: polroots64s
What it does:
  Finds the distinct roots in GF(p) of a polynomial by extracting zero roots, intersecting with x^(p-1)-1, and recursively splitting the linear-factor part.
Inputs:
  - f: degree-d coefficient array, modified in place.
  - d: degree.
  - R: root output array.
  - W: scratch workspace.
  - p: modulus.
Outputs:
  - Returns the number of roots written to R.
Example:
  int nroots = polroots64s(f,d,roots,W,101);
*/
int polroots64s( LONG * f, int d, LONG * R, LONG *W, LONG p )
{
   int i, da, dg; LONG A[2]; extern ULONG seed,mult;
   //printf("roots: deg(f)=%d\n",d);
    // printf("f := "); polprint64s(f,d);
   for( i=0; i<d && f[i]==0; i++ );
   if( i>0 ) { R[0]=0; return( 1 + polroots64s(f+i,d-i,R+1,W,p) ); }
   if( f[d]!=1 ) monic64s(f,d,p);
   A[1] = 1;
   A[0] = 0;
   da = polpowmod64s( A, p-1, f, 1, d, W, W+d, p );    // W = x^(p-1) mod f
   //printf("da = %d, a := ",da); polprint64s(W,da);
   if( da==0 && W[0]==1 ) dg = d; // f is all linear factors
   else { W[0] = sub64s(W[0],1,p); dg = polgcd64s( f, W, d, da, p ); }   // f = gcd(f,W-1)
   //printf("g := "); polprint64s(f,dg);
   if( dg==0 ) return 0;
   seed = 1;
   mult = 6364136223846793003ll;
   polsplit64s( f, dg, R, W, p );
   return dg; // number of roots in R
}

/*
ROOT FINDING WRAPPER FOR MAPLE  (Monagan's polroots64s).

f is given by its d+1 coefficients in ASCENDING order, f[k] the coefficient of
x^k, in any integer representatives (they are reduced into [0,p) here).
The roots in GF(p) are written to
rootsOut in increasing order and their count to info[0].  The caller's f is NOT
modified: polroots64s monics its input, writes the gcd back over it and
recurses on f+i for the zero root, so the coefficients are copied to scratch
first.

WORKSPACE.  polsplit64s documents "a scratch array of size at least 3*d", and
polroots64s's own polpowmod64s call needs C of size d plus R of size 2*db-1,
which is the same 3*d.  The default below is 3*d+8.  wsize overrides it; pass
0 for the default.

Return value:  0  success
              <0  bad arguments
               1  outLen too small; info[0] holds the number of roots found
*/

/*
Function: cppPolRoots
What it does:
  C-compatible wrapper around polroots64s that preserves the caller's polynomial, normalizes coefficients, sizes workspace, and sorts the roots for Maple.
Inputs:
  - d,f: advertised degree and coefficient array.
  - p: prime modulus.
  - wsize: requested scratch size, or 0 for the default.
  - rootsOut/outLen: root output buffer.
  - info/infoLen: metadata; info[0] receives root count.
Outputs:
  - Returns 0 on success, 1 if rootsOut is too small, or a negative argument/internal error code.
Example:
  LONG info[1]; int rc=cppPolRoots(d,f,101,0,d,roots,d?1:1,info);
*/
extern "C" int cppPolRoots(int d,
    const LONG *f,
    const LONG p,
    int wsize,
    int outLen,
    LONG *rootsOut,
    int infoLen,
    LONG *info)
{

// INITIAL CHECKS:

if (!f || !rootsOut || !info) {
    return -1;
}
if (infoLen < 1 || outLen < 0) {
    return -2;
}
if (p < 3) {
    return -3;                      // polsplit64s uses (p-1)/2
}

info[0] = 0;
if (d < 1) {
    return 0;
}

// SCRATCH COPY OF f, REDUCED INTO [0,p) AND TRIMMED TO ITS TRUE DEGREE

std::vector<LONG> F(f, f+d+1);
for (int i = 0; i <= d; ++i) {
    LONG v = F[i] % p;
    F[i] = v + ((v>>63)&p);
}
int df = d;
while (df > 0 && F[df] == 0) df--;
if (df < 1) {
    return 0;                       // f is a nonzero constant
}

// f = c*x^k has 0 as its only root.  polroots64s would deflate down to a
// constant and then divide by it, so take that case here.

int lo = 0;
while (lo < df && F[lo] == 0) lo++;
if (lo == df) {
    info[0] = 1;
    if (outLen < 1) return 1;
    rootsOut[0] = 0;
    return 0;
}

// WORKSPACE

int need = 3*df + 8;
if (wsize < need) wsize = need;
std::vector<LONG> W((size_t)wsize, 0);
std::vector<LONG> R((size_t)df+1, 0);

int n = polroots64s(F.data(), df, R.data(), W.data(), p);
if (n < 0 || n > df) {
    return -4;
}

info[0] = n;
if (n > outLen) {
    return 1;
}

std::sort(R.begin(), R.begin()+n);   // match the ordering of Roots(f) mod p
for (int i = 0; i < n; ++i) {
    rootsOut[i] = R[i];
}

return 0;
}
