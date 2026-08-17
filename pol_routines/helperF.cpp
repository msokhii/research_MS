#include<iostream>
#include<vector>
#include<cstdint>
#include"integerMath.h"

using namespace std; 

using LONG=int_fast64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

void makeMonic(vector<LONG> &v,int deg,LONG p){
	if(deg<0) return;
	LONG LC=v[deg];
	if(LC==1) return;
	LONG invLC=modinv64b(LC,p);
	for(int i=0;i<=deg;i++){
		v[i]=mul64b(v[i],invLC,p);
	}
}

vector<LONG> slicePoly(const vector<LONG>& v, int deg, LONG p){
    if(deg<0){return {};}            // represent 0 as empty
    vector<LONG>out(deg+1,0);
    for(int i=0;i<=deg;i++){
        LONG c=v[i]%p;
        if(c<0){c+=p;}
        out[i]=c;
    }
    // Trim trailing zeroes.
    while(!out.empty() && out.back()==0){out.pop_back();}
    return out;
}
