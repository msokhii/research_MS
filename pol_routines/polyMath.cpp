// Compile with g++ -I .Iheader_FILES pol_ALGO/polyMath.cpp pol_ALGO/integerMath.cpp

#include<iostream>
#include"integerMath.h"
#include<vector> 
#include<unordered_map>
#include<cstdint>
#include"algorithm"
#include"int128g.hpp"

using namespace std;
// extern long long GLOBALCPUMUL;

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

/*
This struct is for pGCDEXTFULL.
*/
struct GCDEX{
	vector<LONG> r;
	vector<LONG> s;
	vector<LONG> t;
	int degR;	
	int degS;
	int degT;
};

/*
This struct is for rational function reconstruction.
*/
struct pairRFR{
	vector<LONG> r;
	vector<LONG> t;
	int degR;
	int degT;
	int flag;
};

/*
This struct returns all values of r,s,t at each iteration.
*/
struct GCDEXHIST{
	GCDEX g;
	vector<vector<LONG>> rTrace;
	vector<vector<LONG>> sTrace;
	vector<vector<LONG>> tTrace;
	vector<int> degRT;
	vector<int> degST;
	vector<int> degTT;
};

/******************************************************************************************/
/* Fast CPU routines                                                                      */
/******************************************************************************************/

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

/******************************************************************************************/
/* Polynomial routines                                                                    */
/******************************************************************************************/

// Prints out the coefficient array of a polynomial from 
// low to high (a0+a1x^1+a2x^2+...+anx^n).

// Computes a^n mod p. Here, 0<=a<p<2^63.

/* 

LONG powmodP64(LONG a,LONG n,LONG p,recint P){
	LONG r;
	LONG s;
	a+=(a>>63)&p;
	if(n==0){return 1;}
	if(n==1){return a;}
	for(r=1,s=a;n>0;n/=2){
		if(n&1){
			r=mulrec64(r,s,P);
			s=mulrec64(s,s,P);
		}
	}
}
	
*/

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

// Returns a pair containing the new vector c=(a*b) mod p
// and the degree of c. 

pair<vector<LONG>,int> pMULNEW64(const vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p){
	vector<LONG> c;
	if(degA<0 || degB<0) return {c,-1};
	int degC=degA+degB;
	c.resize(degC+1,0);
	for(int i=0;i<=degA;i++){
		for(int j=0;j<=degB;j++){
			LONG prod=mul64b(a[i],b[j],p);
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
			temp[i]=mul64b(a[i],x,p);
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
			a[i]=mul64b(a[i],x,p);
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
	t=mul64b(bVal,b[0],p);
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
	t=mul64b(aVal,b[degB],p);
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
        r=add64b(mul64b(r,alpha,p),a[k],p);
    }
    return r;
}

LONG pEVAL64(LONG *a,int d,LONG x,const LONG p){
	int i;
	LONG r;
	if(d==-1){return 0;}
	for(r=a[d],i=d-1;i>=0;i--){
		r=add64b(a[i],mul64b(x,r,p),p);
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
    	for(int i=0;i<=degA;i++) a[i] = mul64b(a[i], inv0, p);
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
		LONG prod1=mul64b(LR,invLTB,p);
		for(int j=0;j<=degB;j++){
			LONG prod2=mul64b(prod1,b[j],p);
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
    for (int i = 0; i <= degA; i++) a[i] = mul64b(a[i], inv0, p);

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
				a[k]=sub64b(a[k],mul64b(t,b[k],p),p);
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
            t=mul64b(t,inv,p);
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
            t=mul64b(t,inv,p);
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
		a[i]=mul64b(invTerm,a[i],p);
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

// This computes the gcd(a(x),b(x)) and puts the GCD i.e. g
// in A and returns degree of g. Both original a and b are 
// updated and-or destroyed. 

/*
int polGCD64(vector<LONG> &a, vector<LONG> &b, int degA, int degB, const LONG p) {
    int degR;
    vector<LONG> c;
	vector<LONG> d;
    LONG u;
	LONG aVal;
	LONG bVal;

    // Division by 0 polynomial.
    if (degB<0){
        cout<<"DIV by 0.\n";
        exit(1);
    }

    // Switches pointers internally but destroys a,b.
    c.swap(a);
    d.swap(b);

    // Make sure degC>=degD.
    if(degA<degB){
        swap(c,d);
        swap(degA,degB);
    }

    while(true){
        // Special case: quotient must be linear when degA=degB+1 (and degB>0).
        if(degB>0 && degA-degB==1){
            u=modinv64b(d[degB],p);
            aVal=mul64b(c[degA],u,p);
            bVal=mul64b(aVal,d[degB-1],p);
            bVal=mul64b(u,sub64b(c[degA-1],bVal,p),p); // quotient=ax+b.
			degR=polSUBMUL64(c,d,aVal,bVal,degA,degB,p); // c=c-(ax+b)d.
			if(degR>=degB){cout << "FAIL.\n";}
        }else{
            // General case: compute remainder of c by d (in-place in c)
            degR=polDIVIP64(c,d,degA,degB,p);
        }

        // The remainder is zero -> gcd is d.
        if (degR<0){
            a.swap(d);              // move gcd into a.
            polMAKEMONIC64(a,p);
            return degB;
        }

        // Else continue the algorithm.
        swap(c,d); 
        degA=degB;
        degB=degR;
    }
}

*/

// Slow version of the extended euclidean alogorithm. Returns 
// all vectors {r,s,t,degrees of all} where r vector contains the gcd of a(x),b(x).
// This is monic.

/* 
GCDEX pGCDEXFULLSLOW(vector<LONG> &r0,vector<LONG> &r1,int degr0,int degr1,const LONG p){
	vector<LONG> s0{1};
	vector<LONG> s1;
	vector<LONG> t0;
	vector<LONG> t1{1};
	int degS0=0;
	int degS1=-1;
	int degT0=-1;
	int degT1=1;
	while(degr1!=-1){
		int degB=degr1;
		pair<int,int> QR=pDIVDEG(r0,r1,degr0,degr1,p);
		int degR=QR.second;
		int degQ=QR.first;
		vector<LONG> q;
		if(degQ>=0){
			q.resize(degQ+1);
			for(int k=0;k<=degQ;k++){
				q[k]=r0[degB+k];
			}
		}
		auto [qs1,dqs1]=pMULNEW64(q,s1,degQ,degS1,p);
        auto [s2,ds2]=pSUBNEW64(polfms64s(LONG *A, LONG *B, LONG *C, int da, int db, int dc, LONG p),qs1,degS0,dqs1,p);
        auto [qt1,dqt1]=pMULNEW64(q,t1,degQ,degT1,p);
        auto [t2,dt2]=pSUBNEW64(t0,qt1,degT0,dqt1,p);
		if(degR==-1){r0.clear();}
		else{r0.resize(degR+1);}
		r0.swap(r1);
		degr0=degB;
		degr1=degR;
		s0.swap(s1); degS0 = degS1;
        s1.swap(s2); degS1 = ds2;
        t0.swap(t1); degT0 = degT1;
        t1.swap(t2); degT1 = dt2;
	}
	polMAKEMONIC64(r0,p);
	return{r0,s0,t0,degr0,degS0,degT0};
}

*/

// Fast version of extended euclidean algorithm. This is done in place.
// Also monic.

/* 
GCDEX pGCDEXFULLFAST(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p){
	if(degA<0 || degB<0){
		cout<<"INPUTS MUST BE NON-ZERO.\n";
		exit(1);
	}
	LONG u;
	LONG aVal;
	LONG bVal;
	int degR;
	int degS;
	int degT;
	int degQ;
	int maxCoeff=max(degA+1,degB+1);
	vector<LONG> s1(maxCoeff,0);
	vector<LONG> s2(maxCoeff,0);
	vector<LONG> t1(maxCoeff,0);
	vector<LONG> t2(maxCoeff,0);
	/*
	s1=1,s2=0
	t1=0,t2=1
	*/
    /*
	s1[0]=1;
	t2[0]=1;
	int degs1=0;
	int degs2=-1;
	int degt1=-1;
	int degt2=0;
	vector<LONG> r1;
	vector<LONG> r2;
	//r1 gets a and r2 gets b. No copying only switching pointers.
	r1.swap(a);
	r2.swap(b);
	// We swap everything if deg(a)<deg(b).
	if(degA<degB){
        swap(r1,r2);
        swap(degA,degB);
        swap(s1,s2);
		swap(degs1,degs2);
        swap(t1,t2); 
		swap(degt1,degt2);
    }
	while(true){
		if(degB>0 && degA-degB==1){
			u=modinv64b(r2[degB],p);
			aVal=mul64b(r1[degA],u,p);
			bVal=mul64b(aVal,r2[degB-1],p);
			bVal=mul64b(u,sub64b(r1[degA-1],bVal,p),p);
			degR=polSUBMUL64(r1,r2,aVal,bVal,degA,degB,p);
			degS=polSUBMUL64(s1,s2,aVal,bVal,degs1,degs2,p);
			degT=polSUBMUL64(t1,t2,aVal,bVal,degt1,degt2,p);
		}
		else{
			degR=polDIVIP64(r1,r2,degA,degB,p);
			degQ=degA-degB;
			vector<LONG> q;
			q.resize(degQ+1,0);
			for(int i=0;i<=degQ;i++){
				q[i]=r1[degB+i];
			}
			vector<LONG> tmpS=s2;
			vector<LONG> tmpT=t2;
			int degtmpS=degs2;
			int degtmpT=degt2;
			degtmpS=pMULIP64(tmpS,q,degtmpS,degQ,p);
			degS=pSUBIP64(s1,tmpS,degs1,degtmpS,p);
			degtmpT=pMULIP64(tmpT,q,degtmpT,degQ,p);
			degT=pSUBIP64(t1,tmpT,degt1,degtmpT,p);
		}
		if(degR<0){
			GCDEX ret;
			ret.r=r2;
			ret.s=s2;
			ret.t=t2;
			ret.degR=degB;
			ret.degS=degs2;
			ret.degT=degt2;
			if(ret.degR>=0){
				LONG LC=ret.r[ret.degR];
				LC%=p;
				if(LC<0){LC+=p;}
				if(LC!=1){
					u=modinv64b(LC,p);
					polSCMULIP64(ret.r,u,ret.degR,p);
					if(ret.degS>=0){
						polSCMULIP64(ret.s,u,ret.degS,p);
					}
					if(ret.degT>=0){
						polSCMULIP64(ret.t,u,ret.degT,p);
					}
				}
			}
			return ret;
		}
		swap(r1,r2);
		degA=degB;
		degB=degR;
		swap(s1,s2);
		int tempS=degs2;
		degs2=degS;
		degs1=tempS;
		swap(t1,t2);
		int tempT=degt2;
		degt2=degT;
		degt1=tempT;
	}
}

*/

// Returns the struct GCDEXHIST that contains all values of 
// r,s and t for each and every iteration.
/* 
GCDEXHIST pGCDEXSTORE64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p){	
	if(degA<0 || degB<0 ){
		cout<<"INPUTS MUST BE NON-ZERO.\n";
		exit(1);
	}
	LONG u;
	LONG aVal;
	LONG bVal;
	int degR;
	int degQ;
	int degS;
	int degT;
	vector<LONG> r1;
	vector<LONG> r2;
	r1=a;
	r2=b;
	int maxCoeff=max(degA+1,degB+1);
	vector<LONG> s1(maxCoeff,0);
	vector<LONG> s2(maxCoeff,0);
	vector<LONG> t1(maxCoeff,0);
	vector<LONG> t2(maxCoeff,0);
	s1[0]=1;
	t2[0]=1;
	int degS1=0;
	int degS2=-1;
	int degT1=-1;
	int degT2=0;
	if(degA<degB){
        swap(r1,r2);
        swap(degA,degB);
        swap(s1,s2);
		swap(degS1,degS2);
        swap(t1,t2); 
		swap(degT1,degT2);
    }
	GCDEXHIST res;
	// [&] is capture by reference.
	// For now, I am not storing the degrees.
	auto triple=[&](const vector<LONG> &rr,const vector<LONG> &ss,const vector<LONG> &tt,int DR,int DS,int DT)
	{
		res.rTrace.push_back(slicePoly(rr,DR,p));
		res.sTrace.push_back(slicePoly(ss,DS,p));
		res.tTrace.push_back(slicePoly(tt,DT,p));
	};
	// These are the first two remainders.
	triple(r1,s1,t1,degA,degS1,degT1);
	triple(r2,s2,t2,degB,degS2,degT2);
	while(true){
		if(degB>0 && degA-degB==1){
			u=modinv64b(r2[degB],p);
			aVal=mul64b(r1[degA],u,p);
			bVal=mul64b(aVal,r2[degB-1],p);
			bVal=mul64b(u,sub64b(r1[degA-1],bVal,p),p);
			degR=polSUBMUL64(r1,r2,aVal,bVal,degA,degB,p);
			degS=polSUBMUL64(s1,s2,aVal,bVal,degS1,degS2,p);
			degT=polSUBMUL64(t1,t2,aVal,bVal,degT1,degT2,p);
		}
		else{
			degR=polDIVIP64(r1,r2,degA,degB,p);
			degQ=degA-degB;
			vector<LONG> q;
			q.resize(degQ+1,0);
			for(int i=0;i<=degQ;i++){
				q[i]=r1[degB+i];
			}
			vector<LONG> tmpS=s2;
			vector<LONG> tmpT=t2;
			int degtmpS=degS2;
			int degtmpT=degT2;
			degtmpS=pMULIP64(tmpS,q,degtmpS,degQ,p);
			degS=pSUBIP64(s1,tmpS,degS1,degtmpS,p);
			degtmpT=pMULIP64(tmpT,q,degtmpT,degQ,p);
			degT=pSUBIP64(t1,tmpT,degT1,degtmpT,p);
		}
		if(degR>=0){
			triple(r1,s1,t1,degA,degS1,degT1);
		}
		if(degR<0){
			GCDEX ret;
			ret.r=r2;
			ret.s=s2;
			ret.t=t2;
			ret.degR=degB;
			ret.degS=degS2;
			ret.degT=degT2;
			if(ret.degR>=0){
				LONG LC=ret.r[ret.degR];
				LC%=p;
				if(LC<0){LC+=p;}
				if(LC!=1){
					u=modinv64b(LC,p);
					polSCMULIP64(ret.r,u,ret.degR,p);
					if(ret.degS>=0){
						polSCMULIP64(ret.s,u,ret.degS,p);
					}
					if(ret.degT>=0){
						polSCMULIP64(ret.t,u,ret.degT,p);
					}
				}
			}
			res.g=ret;
			triple(res.g.r,res.g.s,res.g.t,res.g.degR,res.g.degS,res.g.degT);
			return res;
		}
		swap(r1,r2);
		degA=degB;
		degB=degR;
		swap(s1,s2);
		int tempS=degS2;
		degS2=degS;
		degS1=tempS;
		swap(t1,t2);
		int tempT=degT2;
		degT2=degT;
		degT1=tempT;
	}
}
*/

static inline void makeDenMonicIP64(vector<LONG> &num, int degNum,
                                    vector<LONG> &den, int degDen,
                                    const LONG p){
    if(degDen < 0) return;

    LONG lc = den[degDen] % p;
    if(lc < 0) lc += p;

    if(lc != 1){
        LONG inv = modinv64b(lc, p);
        polSCMULIP64(den, inv, degDen, p);
        if(degNum >= 0){
            polSCMULIP64(num, inv, degNum, p);
        }
    }
}

static inline bool exactDivIP64(vector<LONG> &num, int &degNum,
                                const vector<LONG> &den, int degDen,
                                const LONG p){
    if(degDen < 0) return false;
    if(degNum < degDen) return false;

    auto QR = pDIVDEG(num, den, degNum, degDen, p);
    int degQ = QR.first;
    int degR = QR.second;

    if(degR != -1){
        return false;
    }

    if(degQ < 0){
        num.clear();
        degNum = -1;
        return true;
    }

    // quotient lives in num[degDen .. degDen+degQ]
    for(int k = 0; k <= degQ; k++){
        num[k] = num[degDen + k];
    }
    num.resize(degQ + 1);
    degNum = degQ;
    return true;
}

/* 
pairRFR ratReconFast(const vector<LONG> &m,
                     const vector<LONG> &u,
                     int degM,
                     int degU,
                     int N,
                     int D,
                     const LONG p){

    if(degM < 0 || degU < 0){
        return {{},{},-1,-1,-10};
    }

    auto boundCheck = [&](int degR, int degT) -> bool {
        if(degR > N) return false;
        if(D < 0) return true;
        return degT <= D;
    };

    // Working copies trimmed once
    vector<LONG> r1(m.begin(), m.begin() + degM + 1);
    vector<LONG> r2(u.begin(), u.begin() + degU + 1);

    int degA = degM;
    int degB = degU;

    int maxCoeff = max(degM + 1, degU + 1);

    vector<LONG> t1(maxCoeff, 0);
    vector<LONG> t2(maxCoeff, 0);
    t2[0] = 1;

    int degT1 = -1;
    int degT2 = 0;

    // Reusable buffers
    vector<LONG> q(maxCoeff, 0);
    vector<LONG> tmpT(maxCoeff, 0);

    if(degA < degB){
        swap(r1, r2);
        swap(degA, degB);
        swap(t1, t2);
        swap(degT1, degT2);
    }

    while(degB != -1){

        if(boundCheck(degB, degT2) && degT2 >= 0){
            vector<LONG> num(r2.begin(), r2.begin() + degB + 1);
            vector<LONG> den(t2.begin(), t2.begin() + degT2 + 1);

            makeDenMonicIP64(num, degB, den, degT2, p);

            pairRFR res;
            res.r = move(num);
            res.t = move(den);
            res.degR = degB;
            res.degT = degT2;
            res.flag = 0;
            return res;
        }

        LONG uInv, aVal, bVal;
        int degR, degQ, degT;

        if(degB > 0 && degA - degB == 1){
            uInv = modinv64b(r2[degB], p);
            aVal = mul64b(r1[degA], uInv, p);
            bVal = mul64b(aVal, r2[degB - 1], p);
            bVal = mul64b(uInv, sub64b(r1[degA - 1], bVal, p), p);

            degR = polSUBMUL64(r1, r2, aVal, bVal, degA, degB, p);
            degT = polSUBMUL64(t1, t2, aVal, bVal, degT1, degT2, p);
        }
        else{
            degR = polDIVIP64(r1, r2, degA, degB, p);
            degQ = degA - degB;

            for(int i = 0; i <= degQ; i++){
                q[i] = r1[degB + i];
            }

            if(degT2 >= 0){
                for(int i = 0; i <= degT2; i++){
                    tmpT[i] = t2[i];
                }

                int degTmpT = degT2;
                degTmpT = pMULIP64(tmpT, q, degTmpT, degQ, p);
                degT = pSUBIP64(t1, tmpT, degT1, degTmpT, p);
            }
            else{
                degT = degT1;
            }
        }

        if(degR < 0){
            break;
        }

        swap(r1, r2);
        degA = degB;
        degB = degR;

        swap(t1, t2);
        int oldDegT2 = degT2;
        degT2 = degT;
        degT1 = oldDegT2;
    }

    return {{},{},-1,-1,-20};
}

/* 
pairRFR ratRecon(const vector<LONG> &m,
                 const vector<LONG> &u,
                 int degM,
                 int degU,
                 int N,
                 int D,
                 const LONG p){

    if(degM < 0 || degU < 0){
        return {{},{},-1,-1,-10};
    }

    auto boundCheck = [&](int degR, int degT) -> bool {
        if(degR > N) return false;
        if(D < 0) return true;
        return degT <= D;
    };

    // Trimmed working copies of m and u
    vector<LONG> r1(m.begin(), m.begin() + degM + 1);
    vector<LONG> r2(u.begin(), u.begin() + degU + 1);
    vector<LONG> mTrim = r1;

    int degA = degM;
    int degB = degU;

    int maxCoeff = max(degM + 1, degU + 1);

    vector<LONG> t1(maxCoeff, 0);
    vector<LONG> t2(maxCoeff, 0);
    t2[0] = 1;

    int degT1 = -1;
    int degT2 = 0;

    // Reusable buffers for Euclid loop
    vector<LONG> q(maxCoeff, 0);
    vector<LONG> tmpT(maxCoeff, 0);

    if(degA < degB){
        swap(r1, r2);
        swap(degA, degB);
        swap(t1, t2);
        swap(degT1, degT2);
    }

    while(degB != -1){

        if(boundCheck(degB, degT2) && degT2 >= 0){

            vector<LONG> num(r2.begin(), r2.begin() + degB + 1);
            vector<LONG> den(t2.begin(), t2.begin() + degT2 + 1);
            int degNum = degB;
            int degDen = degT2;

            // Reduce gcd(num, den) if needed, using ordinary gcd only
            {
                vector<LONG> g1 = num;
                vector<LONG> g2 = den;
                int degG = polGCD64(g1, g2, degNum, degDen, p);

                if(degG > 0){
                    g1.resize(degG + 1);

                    if(!exactDivIP64(num, degNum, g1, degG, p)){
                        return {{},{},-1,-1,-31};
                    }
                    if(!exactDivIP64(den, degDen, g1, degG, p)){
                        return {{},{},-1,-1,-31};
                    }
                }
            }

            // Require gcd(den, m) = 1, again using ordinary gcd only
            {
                vector<LONG> dcopy = den;
                vector<LONG> mcopy = mTrim;
                int degGdm = polGCD64(dcopy, mcopy, degDen, degM, p);

                if(degGdm <= 0){
                    makeDenMonicIP64(num, degNum, den, degDen, p);

                    pairRFR res;
                    res.r = move(num);
                    res.t = move(den);
                    res.degR = degNum;
                    res.degT = degDen;
                    res.flag = 0;
                    return res;
                }
            }
        }

        LONG uInv, aVal, bVal;
        int degR, degQ, degT;

        if(degB > 0 && degA - degB == 1){
            uInv = modinv64b(r2[degB], p);
            aVal = mul64b(r1[degA], uInv, p);
            bVal = mul64b(aVal, r2[degB - 1], p);
            bVal = mul64b(uInv, sub64b(r1[degA - 1], bVal, p), p);

            degR = polSUBMUL64(r1, r2, aVal, bVal, degA, degB, p);
            degT = polSUBMUL64(t1, t2, aVal, bVal, degT1, degT2, p);
        }
        else{
            degR = polDIVIP64(r1, r2, degA, degB, p);
            degQ = degA - degB;

            for(int i = 0; i <= degQ; i++){
                q[i] = r1[degB + i];
            }

            if(degT2 >= 0){
                for(int i = 0; i <= degT2; i++){
                    tmpT[i] = t2[i];
                }

                int degTmpT = degT2;
                degTmpT = pMULIP64(tmpT, q, degTmpT, degQ, p);
                degT = pSUBIP64(t1, tmpT, degT1, degTmpT, p);
            }
            else{
                degT = degT1;
            }
        }

        if(degR < 0){
            break;
        }

        swap(r1, r2);
        degA = degB;
        degB = degR;

        swap(t1, t2);
        int oldDegT2 = degT2;
        degT2 = degT;
        degT1 = oldDegT2;
    }

    return {{},{},-1,-1,-20};
}



/*
pairRFR ratRecon(const vector<LONG> &m,const vector<LONG> &u,int degM,int degU,int N,int D,const LONG p){
    if(degM<0 || degU<0){
        return {{},{},-1,-1,-10};
    }

    // Helper functions (Add them in helper file.)

    // gcd computed using fast extended gcd (returns monic gcd in ret.r)
    auto gcdPoly = [&](const vector<LONG> &A, int degA,const vector<LONG> &B, int degB) -> pair<vector<LONG>,int>
    {
        if(degA < 0 || degB < 0) return {{}, -1};
        vector<LONG> a=A, b=B;
        GCDEX g = pGCDEXFULLFAST(a,b,degA,degB,p);
        return {slicePoly(g.r,g.degR, p),g.degR}; // This is monic. Returns g and deg(g).
    };

    // Exact quotient Q = Num/Den (must divide exactly) otherwise skip, using pDIVDEG.
    auto exactQuotient = [&](const vector<LONG> &Num, int degNum,const vector<LONG> &Den, int degDen) -> pair<vector<LONG>,int>
    {
        if(degDen < 0) return {{}, -1};
        vector<LONG> tmp = Num;
        auto QR = pDIVDEG(tmp, Den, degNum, degDen, p);
        int degQ = QR.first;
        int degR = QR.second;
        if(degR != -1) return {{}, -1}; // not divisible
        vector<LONG> Q;
        if(degQ >= 0){
            Q.resize(degQ+1);
            for(int k=0;k<=degQ;k++){
                Q[k] = tmp[degDen + k];
            }
        }
        return {Q, degQ}; // Extracting quotient coefficients into a vector Q.
    };

    // Make Denom monic (scale both by inv(LC(den)))
    auto makeDenMonic = [&](vector<LONG> &num, int degNum,vector<LONG> &den, int degDen)
    {
        if(degDen < 0) return;
        LONG lc = den[degDen] % p; if(lc < 0) lc += p;
        if(lc != 1){
            LONG inv = modinv64b(lc, p);
            polSCMULIP64(den, inv, degDen, p);
            if(degNum >= 0) polSCMULIP64(num, inv, degNum, p);
        }
    };

    // Euclid Algo.
    vector<LONG> r1=m;
    vector<LONG> r2=u;
    int degA=degM;
    int degB=degU;
    int maxCoeff=max(degA+1,degB+1);
    vector<LONG> t1(maxCoeff,0);
    vector<LONG> t2(maxCoeff,0);
    t2[0]=1;
    int degT1=-1;
    int degT2=0;
    if(degA<degB){
        swap(r1,r2);
        swap(degA,degB);
        swap(t1,t2);
        swap(degT1,degT2);
    }
	We require deg(r)<=N.
	If D>=0 then we also require deg(t)<=D.
	If D<0 then the denominator degree is unbounded.
    auto boundCheck=[&](int degR,int degT){
        if(degR<=N){
            if(D<0) return true;
            return degT<=D; 
        }
        return false;
    };

    // Pre-slice m once for gcd(d,m) checks
    vector<LONG> mTrim = slicePoly(m, degM, p);
    while(degB!=-1){
        if(boundCheck(degB,degT2) && degT2 >= 0){
            vector<LONG> num = slicePoly(r2, degB, p);
            vector<LONG> den = slicePoly(t2, degT2, p);
            int degNum = degB;
            int degDen = degT2;

            // Enforce gcd(num,den)=1 (reduce by gcd if needed)
            {
                auto [g, degG] = gcdPoly(num, degNum, den, degDen);
                if(degG > 0){
                    auto [qn, dqn] = exactQuotient(num, degNum, g, degG);
                    auto [qd, dqd] = exactQuotient(den, degDen, g, degG);
                    if(dqn < 0 || dqd < 0){
                        // should not happen if gcd is correct, treat as fail
                        return {{},{},-1,-1,-31};
                    }
                    num.swap(qn); degNum = dqn;
                    den.swap(qd); degDen = dqd;
                }
                // if degG==0, gcd is 1 already (since gcdPoly returns monic)
            }

            // Enforce gcd(den, m)=1 (den invertible mod m)
            {
                auto [gdm, degGdm] = gcdPoly(den, degDen, mTrim, degM);
                if(degGdm > 0){
                    // reject this candidate and keep searching
                    // (a later remainder might have denom coprime to m)
                    // continue loop:
                    goto CONTINUE_EUCLID;
                }
            }

            // Canonical form i.e. make denominator monic.
            makeDenMonic(num, degNum, den, degDen);
            pairRFR res;
			// No extra copying.
            res.r=move(num);
            res.t=move(den);
            res.degR=degNum;
            res.degT=degDen;
            res.flag=0;
            return res;
        }

		// If something fails above, program continues execution from here.
		// Very similar to monagans implementation.
		CONTINUE_EUCLID:
        	LONG uInv, aVal, bVal;
        	int degR, degQ, degT;

        	if(degB>0 && degA-degB==1){
            	uInv=modinv64b(r2[degB],p);
            	aVal=mul64b(r1[degA],uInv,p);
            	bVal=mul64b(aVal,r2[degB-1],p);
            	bVal=mul64b(uInv,sub64b(r1[degA-1],bVal,p),p);
            	degR=polSUBMUL64(r1,r2,aVal,bVal,degA,degB,p);
            	degT=polSUBMUL64(t1,t2,aVal,bVal,degT1,degT2,p);
        	} else{
            	degR=polDIVIP64(r1,r2,degA,degB,p);
            	degQ=degA-degB;
            	vector<LONG> q(degQ+1,0);
            	for(int i=0;i<=degQ;i++){
                	q[i]=r1[degB+i];
            	}
            	vector<LONG> tmpT=t2;
            	int degTmpT=degT2;
            	degTmpT=pMULIP64(tmpT,q,degTmpT,degQ,p);
            	degT=pSUBIP64(t1,tmpT,degT1,degTmpT,p);
        	}
        	if(degR<0){break;}
        	swap(r1,r2);
        	degA=degB;
        	degB=degR;
			swap(t1,t2);
        	int tmpDegT=degT2;
        	degT2=degT;
        	degT1=tmpDegT;
    	}
    return{{},{},-1,-1,-20};
}

*/

/* 
SHORT SUMMARY OF THE ABOVE ROUTINE: 

1. We run extended euclid on (m,u) which are inputs provided by the user.
2. We track the coefficient t such that r congruent t (mod m).
3. Suppose, we encounter a remainder r and r and t are both small enough:
   a: Reduce (r,t) to coprime.
   b: Require gcd(t,m)=1. 
   c: Make the denominator monic for uniqueness.
   d: Return the final answer r/t.
*/

static inline void makeDenMonicOut64(LONG *num, int degNum,
                                     LONG *den, int degDen,
                                     const LONG p,recint P){
    if(den[degDen]!=1){
        LONG inv=modinv64b(den[degDen],p);
        for(int i=0;i<=degDen;i++){
            den[i]=mulrec64(den[i],inv,P);
        }
        for(int i=0;i<=degNum;i++){
            num[i]=mulrec64(num[i],inv,P);
        }
    }
}

/*
int ratReconFastKernelWS(const vector<LONG> &m,
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
                         int &degTOut,
                         recint P){
    auto boundCheck = [&](int degR, int degT)->bool{
        if(degR > N) return false;
        if(D < 0) return true;
        return degT <= D;
    };

    /* copy inputs into workspace
    for(int i=0;i<=degM;i++){
        W.r1[i] = m[i];
    }
    for(int i=0;i<=degU;i++){
        W.r2[i] = u[i];
    }

    // Test Version.
    std::copy_n(m.data(),degM+1,W.r1.data());
    std::copy_n(u.data(),degU+1,W.r2.data());   
    W.t2[0]=1;

    int degA = degM;
    int degB = degU;
    int degT1 = -1;
    int degT2 = 0;

    if(degA < degB){
        swap(W.r1, W.r2);
        swap(degA, degB);
        swap(W.t1, W.t2);
        swap(degT1, degT2);
    }

    while(degB != -1){

        if(boundCheck(degB, degT2) && degT2 >= 0){
            degROut = degB;
            degTOut = degT2;

            
            for(int i=0;i<=degROut;i++){
                rOut[i] = W.r2[i];
            }
            for(int i=0;i<=degTOut;i++){
                tOut[i] = W.t2[i];
            }
            

            // std::copy_n(W.r2.data(), degROut + 1, rOut);
            // std::copy_n(W.t2.data(), degTOut + 1, tOut);

            // makeDenMonicOut64(rOut,degROut,tOut,degTOut,p,P);
            return 0;
        }

        LONG uInv, aVal, bVal;
        int degR, degQ, degT;

        if(degB > 0 && degA - degB == 1){
            uInv=modinv64b(W.r2[degB], p);
            aVal=mulrec64(W.r1[degA], uInv, P);
       
            bVal=mulrec64(aVal, W.r2[degB-1], P);
            bVal=mulrec64(uInv, sub64b(W.r1[degA-1], bVal, p), P);
            degR=polSUBMUL64P(W.r1.data(),W.r2.data(),aVal,bVal,degA,degB,p,P);
            degT=polSUBMUL64P(W.t1.data(),W.t2.data(),aVal,bVal,degT1,degT2,p,P);
        }
        else{
            degR=polDIVP(W.r1.data(),W.r2.data(),degA,degB,p,P);
            degQ=degA-degB;
            for(int i=0;i<=degQ;i++){
                W.q[i]=W.r1[degB+i];
            }
            
            // std::copy_n(W.r1.data() + degB, degQ + 1, W.q.data());

            if(degT2 >= 0){
                
                for(int i=0;i<=degT2;i++){
                    W.tmpT[i] = W.t2[i];
                }
                

                // std::copy_n(W.t2.data(), degT2 + 1, W.tmpT.data());

                int degTmpT=degT2;
                degTmpT=polMUL64P(W.tmpT.data(),W.q.data(),degTmpT,degQ,p,P);
                degT=pSUBIP64(W.t1.data(),W.tmpT.data(),degT1,degTmpT,p);
                
                // degT = polfms64s(W.t2.data(), W.q.data(), W.t1.data(), degT2, degQ, degT1, p);
            }
            else{
                degT = degT1;
            }
        }

        if(degR < 0){
            break;
        }

        swap(W.r1, W.r2);
        degA = degB;
        degB = degR;

        swap(W.t1, W.t2);
        int oldDegT2 = degT2;
        degT2 = degT;
        degT1 = oldDegT2;
    }

    degROut = -1;
    degTOut = -1;
    return -20;
} 
*/

int ratReconFastKernelWS(const vector<LONG> &m,
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
    int &degTOut,
    recint P){

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
                rOut[i] = mulrec64(rOut[i], lcInv, P);
            }
            for(int i = 0; i <= degTOut; i++){
                tOut[i] = mulrec64(tOut[i], lcInv, P);
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

aVal = mulrec64(W.r1[degA], uInv, P);

bVal = mulrec64(aVal, W.r2[degB - 1], P);
bVal = mulrec64(uInv, sub64b(W.r1[degA - 1], bVal, p), P);

degR = polSUBMUL64P(W.r1.data(), W.r2.data(),
           aVal, bVal, degA, degB, p, P);

degT = polSUBMUL64P(W.t1.data(), W.t2.data(),
           aVal, bVal, degT1, degT2, p, P);
}
else{
// Divide r1 by r2:
// quotient goes into high part of W.r1, remainder stays in low part
degR = polDIVP(W.r1.data(), W.r2.data(), degA, degB, p, P);
degQ = degA - degB;

for(int i = 0; i <= degQ; i++){
W.q[i] = W.r1[degB + i];
}

if(degT2 >= 0){
for(int i = 0; i <= degT2; i++){
W.tmpT[i] = W.t2[i];
}

int degTmpT = degT2;
degTmpT = polMUL64P(W.tmpT.data(), W.q.data(), degTmpT, degQ, p, P);
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
    std::copy_n(m.data(),degM+1,W.r1.data());
    std::copy_n(u.data(),degU+1,W.r2.data());
    W.t2[0]=1;
    int degA=degM;
    int degB=degU;
    int degT1=-1;
    int degT2=0;
    if(degA<degB){
        std::swap(W.r1,W.r2);
        std::swap(degA,degB);
        std::swap(W.t1,W.t2);
        std::swap(degT1,degT2);
    }
    while(degB!=-1){
        if(degB==N){
            degROut=degB;
            degTOut=degT2;
            std::copy_n(W.r2.data(),degROut+1,rOut);
            std::copy_n(W.t2.data(),degTOut+1,tOut);
            return 0;
        }
        LONG uInv,aVal,bVal;
        int degR,degQ,degT;
        if(degB>0 && degA-degB==1){
            uInv=modinv64b(W.r2[degB],p);
            aVal=mul64b(W.r1[degA],uInv,p);
            bVal=mul64b(aVal,W.r2[degB-1],p);
            bVal=mul64b(uInv,sub64b(W.r1[degA-1],bVal,p),p);
            degR=polSUBMUL64(W.r1.data(),W.r2.data(),
                            aVal,bVal,degA,degB,p);
            degT=polSUBMUL64(W.t1.data(),W.t2.data(),
                            aVal,bVal,degT1,degT2,p);
        }
        else{
            degR=polDIVIP64(W.r1.data(),W.r2.data(),degA,degB,p);
            degQ=degA-degB;
            for(int i=0;i<=degQ;i++){
                W.q[i]=W.r1[degB+i];
            }
            if(degT2>=0){
                for(int i=0;i<=degT2;i++){
                    W.tmpT[i]=W.t2[i];
                }

            int degTmpT=degT2;
            degTmpT=pMULIP64(W.tmpT.data(),W.q.data(),degTmpT,degQ,p);
            degT=pSUBIP64(W.t1.data(),W.tmpT.data(),degT1,degTmpT,p);
            }
        else{
            degT=degT1;
        }
        }
        if(degR<0){
            break;
        }
        std::swap(W.r1,W.r2);
        degA=degB;
        degB=degR;
        std::swap(W.t1,W.t2);
        int oldDegT2=degT2;
        degT2=degT;
        degT1=oldDegT2;
    }
    degROut=-1;
    degTOut=-1;
    return -20;
};

int ratRecon2(const vector<LONG> &m,
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
                   int &degTOut,
                   recint P){
    std::copy_n(m.data(),degM+1,W.r1.data());
    std::copy_n(u.data(),degU+1,W.r2.data());
    W.t2[0]=1;
    int degA=degM;
    int degB=degU;
    int degT1=-1;
    int degT2=0;
    if(degA<degB){
        std::swap(W.r1,W.r2);
        std::swap(degA,degB);
        std::swap(W.t1,W.t2);
        std::swap(degT1,degT2);
    }
    while(degB!=-1){
        if(degB==N){
            degROut=degB;
            degTOut=degT2;
            std::copy_n(W.r2.data(),degROut+1,rOut);
            std::copy_n(W.t2.data(),degTOut+1,tOut);
            return 0;
        }
        LONG uInv,aVal,bVal;
        int degR,degQ,degT;
        if(degB>0 && degA-degB==1){
            uInv=modinv64b(W.r2[degB],p);
            aVal=mulrec64(W.r1[degA],uInv,P);
            bVal=mulrec64(aVal,W.r2[degB-1],P);
            bVal=mulrec64(uInv,sub64b(W.r1[degA-1],bVal,p),P);
            degR=polSUBMUL64(W.r1.data(),W.r2.data(),
                            aVal,bVal,degA,degB,p);
            degT=polSUBMUL64(W.t1.data(),W.t2.data(),
                            aVal,bVal,degT1,degT2,p);
        }
        else{
            degR=polDIVIP64(W.r1.data(),W.r2.data(),degA,degB,p);
            degQ=degA-degB;
            for(int i=0;i<=degQ;i++){
                W.q[i]=W.r1[degB+i];
            }
            if(degT2>=0){
                for(int i=0;i<=degT2;i++){
                    W.tmpT[i]=W.t2[i];
                }

            int degTmpT=degT2;
            degTmpT=pMULIP64(W.tmpT.data(),W.q.data(),degTmpT,degQ,p);
            degT=pSUBIP64(W.t1.data(),W.tmpT.data(),degT1,degTmpT,p);
            }
        else{
            degT=degT1;
        }
        }
        if(degR<0){
            break;
        }
        std::swap(W.r1,W.r2);
        degA=degB;
        degB=degR;
        std::swap(W.t1,W.t2);
        int oldDegT2=degT2;
        degT2=degT;
        degT1=oldDegT2;
    }
    degROut=-1;
    degTOut=-1;
    return -20;
};
