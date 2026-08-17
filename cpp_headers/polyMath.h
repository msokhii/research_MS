#pragma once

#include"integerMath.h"
#include"helperF.h"
#include"int128g.hpp"
#include<vector>
#include<cstdint> 
#include<unordered_map>
#include<iostream>

using namespace std;

using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

using INT32  = int32_t;
using INT64  = int64_t;
using UINT32 = uint32_t;
using UINT64 = uint64_t;

struct RatReconFastWS{
    vector<LONG> r1;
    vector<LONG> r2;
    vector<LONG> t1;
    vector<LONG> t2;
    vector<LONG> q;
    vector<LONG> tmpT;

    RatReconFastWS(int degM){
        int n=degM+1;
        r1.resize(n,0);
        r2.resize(n,0);
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

struct GCDEXHIST{
	GCDEX g;
	vector<vector<LONG>> rTrace;
	vector<vector<LONG>> sTrace;
	vector<vector<LONG>> tTrace;
	vector<int> degRT;
	vector<int> degST;
	vector<int> degTT;
};

struct pairRFR{
	vector<LONG> r;
	vector<LONG> t;
	int degR;
	int degT;
	int flag;
};

static inline vector<LONG> genVEC64(const int deg,const LONG p){
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
static inline vector<LONG> vecCOPY64(const vector<LONG> &v){
    vector<LONG> temp; 
	temp=v;
	return temp;
};
static inline void dispVEC64(const vector<LONG> &v){
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
unordered_map<LONG,int> checkPOL64(const vector<LONG> &v,const LONG p);
LONG powmodP64(LONG a,LONG n,LONG p,recint P);
pair<vector<LONG>,int> pADDNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p);
int pADDIP64(vector<LONG> &a,const vector<LONG> &b,int &degA,const int degB,const LONG p);
pair<vector<LONG>,int> pSUBNEW64(const vector<LONG> &a,const vector<LONG> &b,const int degA,const int degB,const LONG p);
int pSUBIP64(LONG *a,
             const LONG *b,
             int degA,
             const int degB,
             const LONG p);
int pMULIP64(LONG *a,
             const LONG* b,
             int degA,
             int degB,
             const LONG p);
int polMUL64P(LONG *a,
              LONG *b,
              int degA,
              int degB,
              LONG p,
              recint P);
int pMULIP64VANDER(vector<LONG> &a, vector<LONG> &b, int degA, int degB, const LONG p);
int polSUBMUL64(LONG *a,
                const LONG *b,
                LONG aVal,
                LONG bVal,
                int degA,
                int degB,
                const LONG p);
int polSUBMUL64P(LONG *a,
                const LONG *b,
                LONG aVal,
                LONG bVal,
                int degA,
                int degB,
                const LONG p,
                recint P);
pair<vector<LONG>,int> pMULNEW64(const vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p);
vector<LONG> polSCMULNEW64(vector<LONG> &a,LONG x,int degA,const LONG p);
void polSCMULIP64(vector<LONG> &a,LONG x,int degA,const LONG p);
int polfms64s(LONG *A, LONG *B, LONG *C, int da, int db, int dc, LONG p);
LONG evalHORN64(vector<LONG> &a,LONG alpha,const LONG p);
LONG pEVAL64(LONG *a,int d,LONG x,const LONG p);
pair<int,int> pDIVDEG(vector<LONG> &a,const vector<LONG> &b,int degA,int degB,const LONG p);
int polDIVIP64(LONG *a,
               const LONG *b,
               int degA,
               int degB,
               const LONG p);
int polDIVP(LONG *a,
            LONG *b,
            int degA,
            int degB,
            LONG p,
            recint P);
void polMAKEMONIC64(vector<LONG> &a,const LONG p);
pair<vector<LONG>,int> polGCDNEW64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
int polGCD64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
GCDEX pGCDEXFULLSLOW(vector<LONG> &r0,vector<LONG> &r1,int degr0,int degr1,const LONG p);
GCDEX pGCDEXFULLFAST(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
GCDEXHIST pGCDEXSTORE64(vector<LONG> &a,vector<LONG> &b,int degA,int degB,const LONG p);
pairRFR ratReconFast(const vector<LONG> &m,
                     const vector<LONG> &u,
                     int degM,
                     int degU,
                     int N,
                     int D,
                     const LONG p);
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
                         recint P);
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
                   int &degTOut);
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
                   recint P);
pairRFR ratRecon(const vector<LONG> &m,const vector<LONG> &u,int degM,int degU,int N,int D,const LONG p);
