#pragma once

#include<cstdint>
#include<vector>
#include"integerMath.h"
#include"polyMath.h"

using namespace std; 

using LONG=int64_t;
using ULNG=uint_fast64_t;
using ULNG128=__uint128_t;

void vSOLVER64(
    const vector<LONG>& m,
    const vector<LONG>& y,
    int n,
    vector<LONG>& a,
    vector<LONG>& M,
    int shift,
    LONG p
);