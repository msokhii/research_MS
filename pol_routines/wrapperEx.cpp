/* 
Compile as follows: 

g++ -O3 -shared -fPIC compile headerFiles -o name of objectFile.so cppFiles etc.
Then run the mapleFile using maple mapleFile.
*/

#include"polyMath.h"
#include"integerMath.h"
#include<algorithm>  
#include<vector> 
#include<cstdint>
using namespace std;

using LONG=int_fast64_t;

/*
lib := "/local-scratch/localhome/mss59/Desktop/research_Works/development/libpoly.so":
pADDIP64 := define_external(
  'pADDIP64_C',
  a_len::integer[4],
  A::ARRAY(0..a_len-1, datatype=integer[8]),
  degA::integer[4],
  b_len::integer[4],
  B::ARRAY(0..b_len-1, datatype=integer[8]),
  degB::integer[4],
  p::integer[8],
  RETURN::integer[4],
  LIB=lib
):

This is an example for a C stype wrapper for the function 
pADDIP64 located in polyMath.cpp. We use the extern keyword 
to disable name mangling as maple expects a C stype ABI.

Maple creates an object A which of type ARRAY which we can then
fill with say arbritrary values which is then passed on to our 
wrapper in terms of a pointer. We then define a cpp container 
vector into which we copy the contents of our pointer that we 
got from maple and then we pass this vector into our original 
cpp function pADDIP64 and return whatever this outputs.
*/

extern "C" int pADDIP64_C(int aLen,LONG *A,int degA,int bLen,const LONG *B,int degB,LONG p){
    if(degA==-1&&degB==-1){return -1;}
    vector<LONG> a,b;
    if(degA>0){a.assign(A,A+(degA+1));}
    if(degB>0){b.assign(B,B+(degB+1));}
    int dA=degA;
    int newDeg=pADDIP64(a,b,dA,degB,p);
    if(newDeg>=aLen){return -999;}

    for(int i=0;i<aLen;i++){
        A[i]=0;
    }
    for(int i=0;i<=newDeg;i++){
        A[i]=a[i];
    }
    return newDeg;
}