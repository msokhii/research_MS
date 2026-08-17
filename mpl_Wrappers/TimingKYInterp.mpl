(* First need to compile run.sh and runWrap2.sh and then call maple _FILE_NAME_ *)

restart:

with(NumberTheory):
with(LinearAlgebra):

libV := "/localhome/mss59/Desktop/research_Works/development/pol_ALGO/fastVSolve.so":
libR := "/localhome/mss59/Desktop/research_Works/development/pol_ALGO/rfr.so":


(* Using 0-based indexing *)
mVSOLVE := define_external(
                           'vSOLVER64C',
                           RR::ARRAY(0..nn-1,datatype=integer[8]),
                           LL::ARRAY(0..nn-1,datatype=integer[8]),
                           nn::integer[4],
                           aa::ARRAY(0..nn-1,datatype=integer[8]),
                           XX::ARRAY(0..nn,datatype=integer[8]),
                           shift::integer[4],
                           pp::integer[8],
                           LIB=libV
):

(* Using 0-based indexing *)
mRATRECON := define_external(
                            'ratRECON_C',
                            mLen::integer[4],
                            degM::integer[4],
                            M::ARRAY(0..mLen-1, datatype=integer[8]),
                            uLen::integer[4],
                            degU::integer[4],
                            U::ARRAY(0..uLen-1, datatype=integer[8]),
                            N::integer[4],
                            DBound::integer[4],
                            p::integer[8],
                            nOLEN::integer[4],
                            nOUT::ARRAY(0..nOLEN-1, datatype=integer[8]),
                            degNOUT::REF(integer[4]),
                            dOLEN::integer[4],
                            dOUT::ARRAY(0..dOLEN-1, datatype=integer[8]),
                            degDOUT::REF(integer[4]),
                            RETURN::integer[4],
                            LIB=libR
);

# Remove type checking
# mRATRECON := subsop(1=(mLen,degM,M,uLen,degU,U,N,DBound,p,nOLEN,nOUT,degNOUT,dOLEN,dOUT,degDOUT),op(mRATRECON));

TIC := proc() option inline;
    return kernelopts(cputime);
end proc:

(* Michaels Code *)
VSolveMap := proc(m, v, p, shift::integer:=0 ) option inline;
local t,i,j,M,x,a,q,r,s;

   t := numelems(v);
   if numelems(m) <> t then
       error "v and m must be the same size";
   fi;

   printf("Maple code Vandermonde solver: t=%d  p=%d\n",t,p);

   M := 1;
   for r in m do
       M := Expand( M*(x-r) ) mod p;
   od;

   a := Vector(t);

   for j to t do
       q := Quo(M,x-m[j],x) mod p;
       r := 1/Eval(q,x=m[j]) mod p;
       s := 0;
       for i to t do
           s := s + v[i]*coeff(q,x,i-1);
       od;
       a[j] := r*s mod p;
       if shift <> 0 then
           r := 1/m[j] mod p;
           r := r &^ shift mod p;
           a[j] := r*a[j] mod p;
       fi;
   od;

   if type(v,list) then
       a := convert(a,list);
   fi;

   return a;
end proc:

VandermondeSolve1 := proc(m, v, p, shift::integer:=0) option inline;
local t,i,R,y,a,M,t0_ext,dt_ext;

    t := numelems(v);
    if numelems(m) <> t then
        error "v and m must be the same size";
    fi;

    R := Array(0..t-1,datatype=integer[8]);
    y := Array(0..t-1,datatype=integer[8]);
    a := Array(0..t-1,datatype=integer[8]);
    M := Array(0..t,datatype=integer[8]);

    for i from 0 to t-1 do
        R[i] := m[i+1];
        y[i] := v[i+1];
        a[i] := 0;
    od:

    for i from 0 to t do
        M[i] := 0;
    od:
    mVSOLVE(R,y,t,a,M,shift,p):

    return [seq(a[i], i=0..t-1)];
end proc:

POLYTOARR0 := proc(poly, var, deg, p) option inline;
local A,i;

    A := Array(0..deg, datatype=integer[8]);
    for i from 0 to deg do
        A[i] := coeff(poly,var,i) mod p;
    od:

    return A;
end proc:

ARRTOPOLYNEW := proc(A,deg,var) option inline;
local i:
    add(A[i]*var^i,i=0..deg);
end proc:

ARRTOPOLY0 := proc(A, deg, var, p)
local i, out;

    if deg < 0 then
        return 0;
    fi;

    out := 0;
    for i from 0 to deg do
        out := out + (A[i] mod p)*var^i;
    od:

    return expand(out) mod p;
end proc:

ArrayDegree0 := proc(A, len) option inline;
local i;
    for i from len-1 by -1 to 0 do
        if A[i] <> 0 then
            return i;
        fi;
    od;
    return -1;
end proc:

lastRatReconRC := 0:

Ratrecon1 := proc(Uin, Min, var,
                  N, DBound, p) option inline;
global lastRatReconRC;
local Upoly, Mpoly, degU, degM, uLen, mLen,
      UArr, MArr, nOLEN, dOLEN, nOUT, dOUT,
      degNOUT, degDOUT, rc, nn, dd, i, t0_ext, dt_ext;

    Upoly := Uin:
    Mpoly := Min:

    degU := degree(Upoly,var):
    degM := degree(Mpoly,var):

    if degU < 0 or degM < 0 then
        lastRatReconRC := -999:
        return FAIL;
    fi;

    uLen := degU+1:
    mLen := degM+1:

    UArr := POLYTOARR0(Upoly,var,degU,p):
    MArr := POLYTOARR0(Mpoly,var,degM,p):

    nOLEN := N+1:
    dOLEN := DBound+1:

    nOUT := Array(0..nOLEN-1, datatype=integer[8]):
    dOUT := Array(0..dOLEN-1, datatype=integer[8]):

    degNOUT := -1:
    degDOUT := -1:

    t0_ext := TIC():
    to 10^6 do:
    rc := mRATRECON(
        mLen, degM, MArr,
        uLen, degU, UArr,
        N, DBound, p,
        nOLEN, nOUT, degNOUT,
        dOLEN, dOUT, degDOUT
    ):
    od:
    dt_ext := TIC()-t0_ext:
    printf("DEGNOUT = %d, DEGDOUT = %d, TIME -> %.9f\n",degNOUT,degDOUT,dt_ext/10^6);
    quit;
    
    lastRatReconRC := rc:

    if rc <> 0 then
        return FAIL;
    fi;

    if degNOUT < 0 or degNOUT > N then
        degNOUT := ArrayDegree0(nOUT, nOLEN);
    fi:

    if degDOUT < 0 or degDOUT > DBound then
        degDOUT := ArrayDegree0(dOUT, dOLEN);
    fi:

    if degNOUT < 0 or degNOUT > N then
        return FAIL;
    fi:

    if degDOUT < 0 or degDOUT > DBound then
        return FAIL;
    fi:

    nn := ARRTOPOLYNEW(nOUT,degNOUT,var):
    dd := ARRTOPOLYNEW(dOUT,degDOUT,var):

    if dd = 0 then
        return FAIL;
    fi:

    return nn/dd;
end proc:

p := prevprime(2^63-1):
# RF := rand():
n := x[1]^5+randpoly([seq(x[i],i=1..2)],terms=10,degree=5) mod p:
d := x[1]^5+randpoly([seq(x[i],i=1..2)],terms=10,degree=5) mod p:

f := n/d:

(* Assuming we know this. *)
termsN := nops(n):
termsD := nops(d):

T := max(termsN,termsD):

(* Assuming we know this. *)
N := degree(n):
DD := degree(d):

test := gcd(n,d):
test;

BLACKBOXF := proc(paramA::list(integer), p)
local sigmaVal,result,i;

    sigmaVal := paramA:
    result := Eval(f,[seq(x[i]=sigmaVal[i],i=1..nops(indets(f)))]) mod p:
    return result:
end proc:

COMPUTEBETA := proc()
local betaVal,i,r;

    r := rand(0..p-1):
    betaVal := table():

    for i from 1 to nops(indets(f))-1 do
        betaVal[i] := r():
    od:

    betaVal := convert(betaVal,list):
    return betaVal:
end proc:

GETT := proc(j,betaList)
local sigmaVal,alphaVal,TVal,i,BLACKBOXT,r,a,seen;

    r := rand(1..p-1):
    TVal := table():

    sigmaVal := [seq(ithprime(i)^j,i=1..nops(indets(f)))]:

    seen := table():
    alphaVal := []:
    while nops(alphaVal) < N+DD+1 do
        a := r():
        if not assigned(seen[a]) then
            seen[a] := true:
            alphaVal := [op(alphaVal), a]:
        fi:
    od:

    BLACKBOXT := proc(paramA::integer)
    local alphaX,result;

        alphaX := paramA:
        result := BLACKBOXF(
            [alphaX,
             seq(betaList[i-1]*alphaX-(betaList[i-1]*sigmaVal[1])+sigmaVal[i],
                 i=2..nops(indets(f)))],
            p
        ):
        return result:
    end proc:

    for i from 1 to (N+DD+1) do
        TVal[i] := BLACKBOXT(alphaVal[i]):
    od:

    TVal := convert(TVal,list):
    return alphaVal, TVal:
end proc:

betaList := COMPUTEBETA():

INTERSTEP := proc(betaList)
global lastRatReconRC;
local alphaVal,TVal,interpVal,ratReconVal,mapRatRecon,M,rr,i,j,
      tries,maxTries,ok,t0,dt;

    print(lastRatReconRC);

    ratReconVal := table():
    mapRatRecon := table():
    maxTries := 10:

    for j from 1 to 2*T do
        ok := false:

        for tries from 1 to maxTries while not ok do
            alphaVal,TVal := GETT(j,betaList):
            interpVal := Interp(alphaVal,TVal,z) mod p:
            M := [seq(z-alphaVal[i],i=1..nops(alphaVal))]:
            M := Expand(convert(M,`*`)) mod p:

            # t0 := TIC():
            rr := Ratrecon1(interpVal,M,z,N,DD,p):
            # dt := TIC()-t0:
            # TADD("rr_cpp_total","rr_cpp_calls", dt):
            # LogCSV("rr_cpp", j, tries, "loop", "Ratrecon1_total", dt):

            # t0 := TIC():
            # mapRatRecon[j] := Ratrecon(interpVal,M,z,N,DD) mod p:
            # dt := TIC()-t0:
            # TADD("rr_map_total","rr_map_calls", dt):
            # LogCSV("rr_map", j, tries, "loop", "Maple_Ratrecon", dt):

            if rr <> FAIL then
                ratReconVal[j] := rr:
                ok := true:
            fi;
        od:

        if not ok then
            error "CPP ratRecon wrapper failed at j=%1 after %2 tries; last rc=%3",
                  j, maxTries, lastRatReconRC;
        fi;
    od:

    ratReconVal := convert(ratReconVal,list):
    mapRatRecon := convert(mapRatRecon,list):

    print(nops(ratReconVal)):
    print(nops(mapRatRecon)):

    return ratReconVal;
end proc:

ratReconVal := INTERSTEP(betaList):

ratNumer := table():
ratDenum := table():

for i from 1 to nops(ratReconVal) do
    ratNumer[i] := numer(ratReconVal[i]):
    ratDenum[i] := denom(ratReconVal[i]):
od:

ratNumer := convert(ratNumer,list):
ratDenum := convert(ratDenum,list):

ROOTTOMON := proc(r::posint)
local mon, i, pr, e, q, rr;

    rr := r:
    mon := 1:

    for i from 1 to nops(indets(f)) do
        pr := ithprime(i):
        e := 0:
        while irem(rr,pr,'q') = 0 do
            rr := q:
            e := e+1:
        od:
        mon := mon*x[i]^e:
    od:

    return mon:
end proc:

ROOTTOEXPS := proc(r::posint)
local exps, i, pr, e, q, rr;

    rr := r:
    exps := Array(1..nops(indets(f))):

    for i from 1 to nops(indets(f)) do
        pr := ithprime(i):
        e := 0:
        while irem(rr,pr,'q') = 0 do
            rr := q:
            e := e+1:
        od:
        exps[i] := e:
    od:

    return [seq(exps[i], i=1..nops(indets(f)))];
end proc:

COEFFFROMROOT := proc(poly, r::posint)
local c, exps, i;

    exps := ROOTTOEXPS(r):
    c := poly:

    for i from 1 to nops(exps) do
        c := coeff(c, x[i], exps[i]):
    od:

    return expand(c) mod p;
end proc:

BUILDSPARSE := proc(coeffs::list, roots::list)
local i, out;

    out := 0:
    for i from 1 to nops(coeffs) do
        out := out + coeffs[i]*ROOTTOMON(roots[i]):
    od:
    return expand(out) mod p:
end proc:

GETV := proc(paramA::list(polynom))
local listEval,result,j;

    result := table():
    listEval := paramA:

    for j from 1 to 2*T do
        result[j] := Eval(listEval[j],z=2^j) mod p:
    od:

    result := convert(result,list):
    return result:
end proc:

vNumer := GETV(ratNumer):
vDenom := GETV(ratDenum):

CREATEHANKELMAT := proc(paramA::list(integer), tdim::posint)
local val,i,j,H;

    val := paramA:
    H := Matrix(tdim,tdim):

    for i from 1 to tdim do
        for j from 1 to tdim do
            H[i,j] := val[i+j-1]:
        od:
    od:

    return H:
end proc:

GETBVEC := proc(paramA::list(integer), tdim::posint)
local values,bVec,i;

    values := paramA:
    bVec := -Vector([seq(values[i],i=tdim+1..2*tdim)]):
    return bVec:
end proc:

GETMONFACTOR := proc(paramA::Matrix, paramB::Vector, tdim::posint)
local HMat,bVec,L,annP,rootP,R,F,i,t0Factor,dtFactor;

    HMat := paramA:
    bVec := paramB:
    R := table():
    F := table():

    L := LinearSolve(HMat,bVec) mod p:

    annP := [z^tdim,seq(z^(i-1)*L[i],i=1..tdim)]:
    annP := add(annP[i],i=1..nops(annP)):

    t0Factor := TIC():
    annP := Factor(annP) mod p:
    dtFactor := TIC()-t0Factor:
    TADD("factor_total","factor_calls", dtFactor):
    LogCSV("factor", "-", "-", "proc", sprintf("tdim=%d", tdim), dtFactor):

    rootP := Roots(annP) mod p:

    for i from 1 to nops(rootP) do
        R[i] := rootP[i][1]:
    od:
    R := convert(R,list):

    for i from 1 to nops(R) do
        F[i] := ifactor(R[i]):
    od:
    F := convert(F,list):

    return R,F:
end proc:

HNumer := CREATEHANKELMAT(vNumer, termsN):
HDenom := CREATEHANKELMAT(vDenom, termsD):

bVecNumer := GETBVEC(vNumer, termsN):
bVecDenom := GETBVEC(vDenom, termsD):

numRoots,monNumFactor := GETMONFACTOR(HNumer,bVecNumer,termsN):
denomRoots,monDenomFactor := GETMONFACTOR(HDenom,bVecDenom,termsD):

numRoots := sort(numRoots):
denomRoots := sort(denomRoots):

numCount := nops(numRoots):
denCount := nops(denomRoots):

CREATEVANMAT := proc(paramA::list,n)
local xVals,i,j,vanMat;

    xVals := paramA:
    vanMat := Matrix(n,n):

    for i from 1 to n do
        for j from 1 to n do
            vanMat[i,j] := (xVals[j]^(i-1)) mod p:
        od:
    od:

    return vanMat:
end proc:

vanMatNumer := CREATEVANMAT(numRoots,numCount):
vanMatDenom := CREATEVANMAT(denomRoots,denCount):

GETBVECVAN := proc(paramA::list(integer),n)
local values,bVec,i;

    values := paramA:
    bVec := Vector([seq(values[i],i=1..n)]):
    return bVec:
end proc:

bVecVanNum := GETBVECVAN(vNumer,numCount):
bVecVanDenom := GETBVECVAN(vDenom,denCount):

NORMALIZEPAIR := proc(numCoeff::list, denCoeff::list, denRoots::list, p::prime)
local i, idx, lam, lamInv, newNum, newDen;

    idx := 0:

    for i from 1 to nops(denRoots) do
        if denRoots[i] = 1 then
            idx := i:
            break:
        fi:
    od:

    if idx = 0 then
        for i from 1 to nops(denCoeff) do
            if denCoeff[i] mod p <> 0 then
                idx := i:
                break:
            fi:
        od:
    fi:

    if idx = 0 then
        error "Cannot normalize -> all denominator coefficients are 0";
    fi:

    lam := denCoeff[idx] mod p:
    if lam = 0 then
        error "Cannot normalize -> all chosen denominator coefficient is 0";
    fi:

    lamInv := 1/lam mod p:

    newNum := [seq(lamInv*numCoeff[i] mod p, i=1..nops(numCoeff))]:
    newDen := [seq(lamInv*denCoeff[i] mod p, i=1..nops(denCoeff))]:

    return newNum, newDen, idx;
end proc:

NORMALIZEORIGINALPAIR := proc(nPoly, dPoly, denRoots::list, idx::posint, p::prime)
local lam, lamInv;

    lam := COEFFFROMROOT(dPoly, denRoots[idx]) mod p:
    if lam = 0 then
        error "Original denominator normalization coefficient is 0";
    fi:

    lamInv := 1/lam mod p:

    return expand(lamInv*nPoly) mod p, expand(lamInv*dPoly) mod p;
end proc:

LVanNum := VandermondeSolve1( [seq(vNumer[i],i=1..numCount)],numRoots, p, 1):

LVanDenom := VandermondeSolve1([seq(vDenom[i],i=1..denCount)],denomRoots, p, 1):




LVanNumMap := VSolveMap([seq(vNumer[i],i=1..numCount)],numRoots, p, 1):


LVanDenomMap := VSolveMap([seq(vDenom[i],i=1..denCount)],denomRoots, p, 1):
LVanNum, LVanDenom, normIdx := NORMALIZEPAIR(LVanNum, LVanDenom, denomRoots, p):

recNum := BUILDSPARSE(LVanNum, numRoots):
recDen := BUILDSPARSE(LVanDenom, denomRoots):

nNorm, dNorm := NORMALIZEORIGINALPAIR(n, d, denomRoots, normIdx, p):

sameRF := evalb(expand(recNum*d - recDen*n) mod p = 0):
sameNum := evalb(expand(recNum - nNorm) mod p = 0):
sameDen := evalb(expand(recDen - dNorm) mod p = 0):

printf("Same rational function? -> %a\n", sameRF):
printf("Recovered numerator matches normalized original? -> %a\n", sameNum):
printf("Recovered denominator matches normalized original? -> %a\n", sameDen):
