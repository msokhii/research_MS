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

inline LONG add64b(LONG a,LONG b,LONG p){
    LONG r=(a+b)-p;
    r+=(r>>63)&p;
    return r;
}

inline LONG sub64b(LONG a,LONG b,LONG p){
    LONG r=(a-b);
    r+=(r>>63)&p;
    return r;
}

inline LONG mul64b(LONG a,LONG b, LONG p){
    ULNG128 res=(ULNG128)a*b;
    ULNG r=(ULNG)(res%p);
    return r;
}

inline LONG neg64s(LONG a,LONG p){ 
    return (a==0)?0:p-a; 
};

// We are assuming 0<=a,b<p for the following routines. 

inline LONG mul64bASM(LONG a,LONG b,LONG p){
    LONG q, r;
    __asm__ __volatile__(           \
    "       mulq    %%rdx           \n\t" \
    "       divq    %4              \n\t" \
    : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
    return r;
}

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

vector<LONG> vecCOPY64(const vector<LONG> &v){
vector<LONG> temp; 
temp=v;
return temp;
};

void vecfill64s( LONG x, LONG *A, int n )
{   int i;
    for( i=0; i<n; i++ ) A[i] = x;
    return;
}

void polcopy64s( LONG *A, int d, LONG *B )
{   int i;
    for( i=0; i<=d; i++) B[i]=A[i];
    return;
}

void monic64s(LONG *A,int d,LONG p) {
    int i; LONG inv;
    if(d<0 || A[d]==1) return;
    inv = modinv64b(A[d],p);
    for(i=0; i<d; i++) A[i] = mul64b(inv,A[i],p);
    A[d] = 1;
    return;
}

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

int polMUL64P(LONG *a,
              LONG *b,
              int degA,
              int degB,
              LONG p,
              recint P){
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
            t=add64b(t,mulrec64(a[i],b[k-i],P),p);
        }
        a[k]=t;
    }
    while(degC>=0 && a[degC]==0){
        degC--;
    }
    return degC;
}

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

int polSUBMUL64P(LONG *a,
                 const LONG *b,
                 LONG aVal,
                 LONG bVal,
                 int degA,
                 int degB,
                 const LONG p,
                 recint P){
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
    t=mulrec64(bVal,b[0],P);
    a[0]=sub64b(a[0],t,p);
    for(i=1;i<=degB;i++){
        t=mulrec64(aVal,b[i-1],P);
        t=add64b(t,mulrec64(bVal,b[i],P),p);
        a[i]=sub64b(a[i],t,p);
    }
    t=mulrec64(aVal,b[degB],P);
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

LONG evalHORN64(vector<LONG>& a,LONG alpha,LONG p){
    LONG r = 0LL;
	for (int k=a.size();k-->0;){
        r=add64b(mul64bASM(r,alpha,p),a[k],p);
    }
    return r;
}

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

int polDIVP(LONG *a,
            LONG *b,
            int degA,
            int degB,
            LONG p,
            recint P){
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
            t=sub64b(t,mulrec64(b[j],a[k-j+degB],P),p);
        }
        if(k>=degB && inv!=1){
            t=mulrec64(t,inv,P);
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

static int polExactQuo(const std::vector<LONG>&A,int degA,const std::vector<LONG>&B,int degB,LONG p,std::vector<LONG>&Q){
    std::vector<LONG> a=A; a.resize(degA+1);
    polDIVIP64(a.data(),B.data(),degA,degB,p);
    int degQ=degA-degB;
    Q.assign(a.begin()+degB,a.begin()+degB+degQ+1);
    while(degQ>=0 && Q[degQ]==0) degQ--;
    Q.resize(degQ<0?1:degQ+1);
    return degQ;
}

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
         //s = mulrec64(u,s,P);
         if( shift!=0 ) {
             u = modinv64b(m[j],p);
             u = powmod64s(u,shift,p);
             //s = mulrec64(u,s,P);
             s = mul64b(u,s,p);
         }
         a[j] = s;
         polmul64s(A,M,M,1,n-1,p); // restore M = M x (x-m[j])
    }
    return;
}

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

int newtonInterpMulRec(LONG* x,
    LONG* y,
    const int n,
    const LONG p,
    recint P){
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
            s=add64b(s,mulrec64(prod,Y[i],P),p);
	        prod=mulrec64(prod,sub64b(X[j],X[i],p),P);
        }
        if(prod==0){
            return -1;
        }
        Y[j]=mulrec64(sub64b(Y[j],s,p),modinv64b(prod,p),P);	
    }
    d=n-1;
    while(d>=0&&y[d]==0){
        d--;
    }
    for(i=1;i<=d;i++){
        for(j=d-i;j<=d-1;j++){
            Y[j]=sub64b(Y[j],mulrec64(X[d-i],Y[j+1],P),p);        
	    }
    }
    return d;
}

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

static std::vector<LONG> ftr_buildM(const LONG *alpha,int n,LONG p){
    std::vector<LONG> M(n+1,0); M[0]=1; int degM=0; LONG lin[2]; lin[1]=1;
    for(int i=0;i<n;i++){ lin[0]=(alpha[i]==0?0:p-alpha[i]); degM=pMULIP64(M.data(),lin,degM,1,p); }
    return M;                                  // degree n
}

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
recint P = recip1(p);

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
*/
 
extern "C" int cppRREF(int n,int m,LONG *B,const LONG p){
    if (!B){
        return -1;
    }
    if (n <= 0 || m <= 0){
        return -2;
    }
    if (p<2){
        return -3;
    }
 
    for (size_t i=0; i<(size_t)n*m;++i){
        LONG v = B[i]%p;
        B[i] = v+((v>>63)&p);
    }
 
    return (int)rref(B,n,m,p);
}