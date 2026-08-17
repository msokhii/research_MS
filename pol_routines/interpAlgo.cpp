#include<vector>
#include<cstdint>
#include<cstdlib>
#include"integerMath.h"
#include"int128g.hpp"

using namespace std;

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
