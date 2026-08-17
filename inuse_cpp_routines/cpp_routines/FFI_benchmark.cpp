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
#include"pol_math.cpp"

using namespace std;
using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

int main(){
    LONG p = 2147483647; //This is 2^31-1.
    recint P = recip1(p);

    int degN = 5;
    int degD = 5;

    /*
    Increase this for a more accurate timing. 
    */
    const int CALLS = 10;
    const int ITER  = 2;
    /*
    home_file_pwd = "/home/msokhi/Desktop/res_MS/rat_recon_mpl/timings/FFI_overhead.txt"
    maple_file_pwd = "/cecm/home/mss59/Desktop/resMaple_MS/rat_recon_mpl/timings/FFI_overhead.txt"
    */
    ofstream logFile("/home/msokhi/Desktop/res_MS/rat_recon_mpl/timings/FFI_overhead.txt");
    if (!logFile) {
        cerr << "Could not find directory.\n";
        return 1;
    }

    logFile << "PRIME -> " << p << "\n";
    logFile << "CALLS -> " << CALLS << "\n";
    logFile << left
            << setw(8)  << "Step"
            << setw(8)  << "degN"
            << setw(8)  << "degD"
            << setw(24) << "NewtonMULREC"
            << setw(24) << "NewtonMUL64"
            << setw(24) << "NewtonFFIcpp"
            << setw(24) << "rrMUL64"
            << setw(24) << "rrFFIcpp"
            << "\n";

    for (int step = 1; step < ITER; ++step) {
        vector<LONG> n(degN + 1,0);
        vector<LONG> d(degD + 1,0);

        for (int i = 0; i <= degN; ++i) {
            LONG temp = rand64s(p);
            while (temp == 0) temp = rand64s(p);
            n[i] = temp;
        }

        for (int j = 0; j <= degD; ++j) {
            LONG temp = rand64s(p);
            while (temp == 0) temp = rand64s(p);
            d[j] = temp;
        }

        if (d[degD] != 1) {
            LONG invTerm = modinv64b(d[degD], p);
            for (int i = 0; i <= degD; ++i) d[i] = mul64b(invTerm, d[i], p);
            for (int i = 0; i <= degN; ++i) n[i] = mul64b(invTerm, n[i], p);
        }

        int m = degN + degD + 1;

        vector<LONG> x(m, 0);
        for (int i = 0; i < m; ++i) x[i] = i + 1;

        vector<LONG> yVals(m, 0);
        for (int i = 0; i < m; ++i) {
            LONG denEval = pEVAL64(d.data(), degD, x[i], p);
            if (denEval == 0) {
                cerr << "Cannot divide by 0.\n";
                return 1;
            }
            LONG numEval = pEVAL64(n.data(), degN, x[i], p);
            yVals[i] = mul64b(numEval, modinv64b(denEval, p), p);
        }
        
        vector<LONG> yWork(m, 0);

        auto cpStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(yVals.begin(), yVals.end(), yWork.begin());
        }
        auto cpStop = chrono::steady_clock::now();
        double copyOnly_us =
        chrono::duration<double, std::micro>(cpStop - cpStart).count() / CALLS;

        int degU_rec = -1;
        auto nRecStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(yVals.begin(), yVals.end(), yWork.begin());
            degU_rec = newtonInterpMulRec(x.data(), yWork.data(), m, p, P);
        }
        auto nRecStop = chrono::steady_clock::now();
        double newtonRecWithCopy_us =
        chrono::duration<double, std::micro>(nRecStop - nRecStart).count() / CALLS;
        double newtonKernelRec_us = newtonRecWithCopy_us - copyOnly_us;

        if (degU_rec < 0) {
            cerr << "newtonInterpMulRec failed.\n";
            return 1;
        }

        copy(yVals.begin(), yVals.end(), yWork.begin());
        degU_rec = newtonInterpMulRec(x.data(), yWork.data(), m, p, P);
        if (degU_rec < 0) {
            cerr << "newtonInterpMulRec failed while building Ucoeff.\n";
            return 1;
        }
        vector<LONG> Ucoeff(yWork.begin(), yWork.begin() + (degU_rec + 1));

        int degU_64 = -1;
        auto n64Start = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(yVals.begin(), yVals.end(), yWork.begin());
            degU_64 = newtonInterpMulNormal(x.data(), yWork.data(), m, p);
        }
        auto n64Stop = chrono::steady_clock::now();
        double newton64WithCopy_us =
        chrono::duration<double, std::micro>(n64Stop - n64Start).count() / CALLS;
        double newtonKernel64_us = newton64WithCopy_us - copyOnly_us;

        vector<LONG> yOutWrap(m, 0);
        int degOutWrap = -1;

        auto nWrapStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            int rc = cppInterp(
                m, x.data(),
                m, yVals.data(),
                p,
                m, yOutWrap.data(),
                &degOutWrap
            );
            if (rc != 0) {
                cerr << "cppInterp failed with rc = " << rc << "\n";
                return 1;
            }
        }
        auto nWrapStop = chrono::steady_clock::now();
        double newtonWrapCPP_us =
        chrono::duration<double, std::micro>(nWrapStop - nWrapStart).count() / CALLS;

        vector<LONG> M(degN + degD + 2, 0);  
        M[0] = 1;
        int degM = mkM(M, x, p); 

        if (degM < 0) {
            cerr << "mkM failed.\n";
            return 1;
        }

        vector<LONG> Mcoeff(M.begin(), M.begin() + (degM + 1));
        vector<LONG> Mwork = Mcoeff;
        vector<LONG> Uwork = Ucoeff;

        auto rrCopyStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(Mcoeff.begin(), Mcoeff.end(), Mwork.begin());
            copy(Ucoeff.begin(), Ucoeff.end(), Uwork.begin());
        }
        auto rrCopyStop = chrono::steady_clock::now();
        double rrCopyOnly_us =
        chrono::duration<double, std::micro>(rrCopyStop - rrCopyStart).count() / CALLS;

        RatReconFastWS W(degM);
        vector<LONG> rOut(degN + 1, 0);
        vector<LONG> tOut(degD + 1, 0);

        int flag = -999;
        int degR = -1;
        int degT = -1;

        auto rrKernelStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            copy(Mcoeff.begin(), Mcoeff.end(), Mwork.begin());
            copy(Ucoeff.begin(), Ucoeff.end(), Uwork.begin());

            degR = -1;
            degT = -1;
            flag = ratReconNormal(
                Mwork,
                Uwork,
                degM,
                degU_rec,
                degN,
                degD,
                p,
                W,
                rOut.data(),
                degR,
                tOut.data(),
                degT
            );

            if (flag != 0) {
                cerr << "ratReconNormal failed with rc = " << flag << "\n";
                return 1;
            }
        }
        auto rrKernelStop = chrono::steady_clock::now();
        double rrKernelWithCopy_us =
        chrono::duration<double, std::micro>(rrKernelStop - rrKernelStart).count() / CALLS;
        double rrKernelFastWS_us = rrKernelWithCopy_us - rrCopyOnly_us;

        vector<LONG> nOutWrap(degN + 1, 0);
        vector<LONG> dOutWrap(degD + 1, 0);
        int degNOutWrap = -1;
        int degDOutWrap = -1;

        auto rrWrapStart = chrono::steady_clock::now();
        for (int i = 0; i < CALLS; ++i) {
            int rc = ratRECON_C(
                degM + 1, degM, Mcoeff.data(),
                degU_rec + 1, degU_rec, Ucoeff.data(),
                degN, degD, p,
                degN + 1, nOutWrap.data(), &degNOutWrap,
                degD + 1, dOutWrap.data(), &degDOutWrap
            );
            if (rc != 0) {
                cerr << "ratRECON_C failed with rc = " << rc << "\n";
                return 1;
            }
        }
        auto rrWrapStop = chrono::steady_clock::now();
        double rrWrapCPP_us =
        chrono::duration<double, std::micro>(rrWrapStop - rrWrapStart).count() / CALLS;
        logFile << left
                << setw(8)  << step
                << setw(8)  << degN
                << setw(8)  << degD
                << setw(24) << newtonKernel64_us
                << setw(24) << newtonKernelRec_us
                << setw(24) << newtonWrapCPP_us
                << setw(24) << rrKernelFastWS_us
                << setw(24) << rrWrapCPP_us
                << "\n";

        degN *= 2;
        degD *= 2;
    }

    logFile.close();
    return 0;
}