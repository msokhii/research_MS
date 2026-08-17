#include"polyMath.h"
#include"integerMath.h"
#include"interpAlgo.h"
#include"rfrWrapper.h"
#include<vector>
#include<iostream>
#include<cstdint>
#include<time.h>
#include<chrono>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std; 

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

long long GLOBALMUL=0;
long long GLOBALMUL64=0;
long long GLOBALCPUMUL=0;
int main(){

    LONG p=4294967291; // This is prevprime(2^32-1) from maple. 
    recint P=recip1(p);
    int degN=5;
    int degD=5;
    const int CALLS=1000; 
    const int ITER=5;

    ofstream logFile("benchMark.txt");
    logFile<<"PRIME -> "<<p<<"\n";
    logFile<<"CALLS -> "<<CALLS<<"\n";
    logFile<<left
        <<setw(10)<<"ITER"
        <<setw(10)<<"degN"
        <<setw(10)<<"degD"
        <<setw(28)<<"avgTimeNewton(mulRec)"
        <<setw(28)<<"avgTimeNewton(mul64)"
        <<setw(28)<<"avgTimeRR(No CPU+mulRec)"
        <<setw(28)<<"avgTimeRR(CPU+mul64)"
        <<setw(28)<<"avgTimeRR(CPU+mulRec)"
        <<setw(22)<<"mulsNewton(mulRec)"
        <<setw(22)<<"mulsNewton(mulSub)"
        <<setw(18)<<"mulsRR"
        << "\n";

    for(int step=1;step<ITER;step++){
        vector<LONG>n(degN+1,0);
        vector<LONG>d(degD+1,0);
        /* 
        Populate the numerator vector.
        */
        for(int i=0;i<degN+1;i++){
            LONG temp=rand64s(p);
            while(temp==0){
                temp=rand64s(p);
            }
            n[i]=temp;
        }
        /* 
        Populate the denominator vector.
        */
        for(int j=0;j<degD+1;j++){
            LONG temp=rand64s(p);
            while(temp==0){
                temp=rand64s(p);
            }
            d[j]=temp;
        }
        /* 
        Make the denominator monic.
        */
        if(d[degD]!=1){
            LONG invTerm;
            invTerm=modinv64b(d[degD],p);
            for(int i=0;i<=degD;i++){
                d[i]=mul64b(invTerm,d[i],p);
            }
            for(int j=0;j<=degN;j++){
                n[j]=mul64b(invTerm,n[j],p);
            }
        }

        vector<LONG>nCopy=n;
        vector<LONG>dCopy=d;
        /* 
        We need degN+degD+1 points to interpolate. Here we make the x vector. 
        */
        int m=degN+degD+1;
        vector<LONG>x(m,0);
        for(int i=0;i<m;i++){
            x[i]=i+1;
        }   
        

        vector<LONG>y(m,0);
        for(int i=0;i<m;i++){
            LONG denEval=pEVAL64(d.data(),degD,x[i],p);
            if(denEval==0){
                return -1;
            }
            LONG numEval=pEVAL64(n.data(),degN,x[i],p);
            y[i]=mul64b(numEval,modinv64b(denEval,p),p);
        }

        vector<LONG>yCopy(m,0);
        copy(y.begin(),y.end(),yCopy.begin());
        int degU=newtonInterpMulRec(x.data(),yCopy.data(),m,p,P);
        
        // Timer for newton interpolation using mulrec64 routine.
        GLOBALMUL=0;
        auto start=chrono::steady_clock::now();
        for(int i=0;i<CALLS;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterpMulRec(x.data(),yCopy.data(),m,p,P);
        };
        auto stop=chrono::steady_clock::now();
        
        // Timer for newton interpolation using mul64b routine.
        GLOBALMUL64=0;
        auto newton2Start=chrono::steady_clock::now();
        for(int i=0;i<CALLS;i++){
            copy(y.begin(),y.end(),yCopy.begin());
            int degU=newtonInterpMulNormal(x.data(),yCopy.data(),m,p);
        }
        auto newton2Stop=chrono::steady_clock::now();
        
        // Timer for copying y into y0 for newton interpolation.
        auto cpStart=chrono::steady_clock::now();
        for(int i=0;i<CALLS;i++){
            copy(y.begin(),y.end(),yCopy.begin());
        }
        auto cpStop=chrono::steady_clock::now();
        double cpTotal=chrono::duration<double,std::micro>(cpStop-cpStart).count();
        double total=chrono::duration<double,std::micro>(stop-start).count();
        double newton2Total=chrono::duration<double,std::micro>(newton2Stop-newton2Start).count();
        double avgTimeCp=cpTotal/CALLS;
        double avgTimeNewton=(total/CALLS)-avgTimeCp;
        double avgTimeNewton2=(newton2Total/CALLS)-avgTimeCp;
        long long newtonMuls=GLOBALMUL/CALLS;
        long long newtonMuls2=GLOBALMUL64/CALLS;
        vector<LONG>M(m+1,0);
        M[0]=1;
        int degM=mkM(M,x,p);

        RatReconFastWS W(degM);
        RatReconFastWS W2(degM);
        RatReconFastWS W3(degM);
        vector<LONG>rOut(m,0);
        vector<LONG>tOut(m,0);
        vector<LONG>rOut2(m,0);
        vector<LONG>tOut2(m,0);
        vector<LONG>rOut3(m,0);
        vector<LONG>tOut3(m,0);
        int degR=-1;
        int degT=-1;
        int flag=-999;
        int degR2=-1;
        int degT2=-1;
        int flag2=-999;
        int degR3=-1;
        int degT3=-1;
        int flag3=-999;
        int degUCP=degU;
        int degNCP=degN;
        int degDCP=degD;
        int mCP=m;
        int degUCP3=degU;
        int degNCP3=degN;
        int degDCP3=degD;
        int mCP3=m;
        vector<LONG> MCP=M;
        vector<LONG> yCP2=y;
        vector<LONG> MCP3=M;
        vector<LONG> yCP3=y;
        GLOBALMUL=0;
        auto start2=chrono::steady_clock::now();
        for(int k=0;k<CALLS;k++){
            flag=ratReconFastKernelWS(M,y,m,degU,
            degN,degD,p,W,rOut.data(),degR,tOut.data(),degT,P);
        }
        auto stop2=chrono::steady_clock::now();
        GLOBALCPUMUL=0;
        GLOBALMUL64=0;
        auto rrNormStart=chrono::steady_clock::now();
        for(int k=0;k<CALLS;k++){
            flag2=ratReconNormal(MCP,yCP2,mCP,degUCP,
            degNCP,degDCP,p,W2,rOut2.data(),degR2,tOut2.data(),degT2);
        }
        auto rrNormStop=chrono::steady_clock::now();
        auto rrNorm2Start=chrono::steady_clock::now();
        for(int k=0;k<CALLS;k++){
            flag3=ratRecon2(MCP3,yCP3,mCP3,degUCP3,
            degNCP3,degDCP3,p,W3,rOut3.data(),degR3,tOut3.data(),degT3,P);
        }
        auto rrNorm2Stop=chrono::steady_clock::now();
        
        long long RRmuls=GLOBALMUL/CALLS;
        long long cpuMULRR=GLOBALCPUMUL/CALLS;
        long long totalRRMuls=RRmuls;
        double total2=chrono::duration<double,std::micro>(stop2-start2).count();
        double rrNormTotal=chrono::duration<double,std::micro>(rrNormStop-rrNormStart).count();
        double rrNorm2Total=chrono::duration<double,std::micro>(rrNorm2Stop-rrNorm2Start).count();
        double avgTimeRR=total2/CALLS;
        double avgTimeRRNorm=rrNormTotal/CALLS;
        double avgTimeRRNorm2=rrNorm2Total/CALLS;
        
        logFile<<left<<
                 setw(10)<<step<<
                 setw(10)<<degN<<
                 setw(10)<<degD<<
                 setw(28)<<avgTimeNewton<<
                 setw(28)<<avgTimeNewton2<<
                 setw(28)<<avgTimeRR<<
                 setw(28)<<avgTimeRRNorm<<
                 setw(28)<<avgTimeRRNorm2<<
                 setw(22)<<newtonMuls<<
                 setw(22)<<newtonMuls2<<
                 setw(18)<<totalRRMuls<<
                 "\n";
                  
        degN*=2;
        degD*=2;
    }
    logFile.close();
    return 0;
}
