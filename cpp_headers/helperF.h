#pragma once

#include"integerMath.h"
#include<vector>
#include<cstdint>

using namespace std;

void makeMonic(vector<LONG> &v,int deg,LONG p);
vector<LONG> slicePoly(const vector<LONG>& v, int deg, LONG p);

