#pragma once

#include <cstdint>
using namespace std;
using LONG = int64_t;

#ifdef __cplusplus
extern "C" {
#endif

int ratRECON_C(int mLen,
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
               int *degDOUT);

#ifdef __cplusplus
}
#endif
