#pragma once

#include"integerMath.h"
#include"helperF.h"
#include"int128g.hpp"

int newtonInterpMulRec(LONG *x,
                       LONG *y,
                       const int n,
                       const LONG p,
                       recint P);
int newtonInterpMulNormal(LONG* x,
                          LONG* y,
                          const int n,
                          const LONG p);
