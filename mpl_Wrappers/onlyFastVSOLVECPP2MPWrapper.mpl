restart:

# ------------------------------------------------------------
# External C Vandermonde solver
# ------------------------------------------------------------
libV := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/wrapOBJ2.so":

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

(* 
There is some confusion with how the wrapper takes in inputs. Fix later.
*)
VandermondeSolve1 := proc(m::{Vector,list}, v::{Vector,list}, p::prime, shift::integer:=0)
local t,i,R,y,a,M;

    t := numelems(v);
    if numelems(m) <> t then
        error "v and m must be the same size";
    fi;

    printf("Quadratic C code Vandermonde solver: t=%d  p=%d\n",t,p);

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

    mVSOLVE(R,y,t,a,M,shift,p);

    return [seq(a[i], i=0..t-1)];
end:

# ------------------------------------------------------------
# Main setup
# ------------------------------------------------------------
with(NumberTheory):
with(LinearAlgebra):

p := prevprime(2^31-1):
printf("Prime chosen: %a\n",p):

n := randpoly([x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]],terms=300) mod p:
d := randpoly([x[1],x[2],x[3],x[4],x[5],x[6]],terms=100) mod p:

printf("NUMERATOR -> \n"):
n;

printf("DENOMINATOR -> \n"):
d;

printf("No. of terms in NUM: %a\nNo. of terms in DENUM: %a\n",nops(n),nops(d)):

f := n/d:
printf("RATIONAL FUNCTION -> \n"):
f;

termsN := nops(n):
termsD := nops(d):

# T is only a bound for how many univariate samples we generate
T := max(termsN,termsD):
printf("Value of T: %a\n",T):

N := degree(n):
DD := degree(d):
printf("Degree of NUM: %a\nDegree of DENUM: %a\n",N,DD):

test := gcd(n,d):
printf("GCD of NUM and DENUM: %a\n",test):

# ------------------------------------------------------------
# Black-box evaluation
# ------------------------------------------------------------
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
local sigmaVal,alphaVal,TVal,i,BLACKBOXT,r;

    r := rand(1..p-1):
    TVal := table():

    sigmaVal := [seq(ithprime(i)^j,i=1..nops(indets(f)))]:
    alphaVal := [seq(r(),i=1..(N+DD+1))]:

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

# ------------------------------------------------------------
# Interpolation + built-in Maple Ratrecon
# ------------------------------------------------------------
INTERSTEP := proc(betaList)
local alphaVal,TVal,interpVal,ratReconVal,M,i,j;

    ratReconVal := table():

    for j from 1 to 2*T do
        alphaVal,TVal := GETT(j,betaList):
        interpVal := Interp(alphaVal,TVal,z) mod p:
        M := [seq(z-alphaVal[i],i=1..nops(alphaVal))]:
        M := Expand(convert(M,`*`)) mod p:
        ratReconVal[j] := Ratrecon(interpVal,M,z,N,DD) mod p:
    od:

    ratReconVal := convert(ratReconVal,list):
    return ratReconVal:
end proc:

ratReconVal := INTERSTEP(betaList):

printf("Rational Functions -> \n"):
ratReconVal;

# ------------------------------------------------------------
# Extract numerators / denominators
# ------------------------------------------------------------
ratNumer := table():
ratDenum := table():

for i from 1 to nops(ratReconVal) do
    ratNumer[i] := numer(ratReconVal[i]) mod p:
    ratDenum[i] := denom(ratReconVal[i]) mod p:
od:

ratNumer := convert(ratNumer,list):
ratDenum := convert(ratDenum,list):

printf("Numerators -> \n"):
ratNumer:

printf("Denominators -> \n"):
ratDenum:

# ------------------------------------------------------------
# Utility: decode root integer to monomial
# ------------------------------------------------------------
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
end:

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
end:

COEFFFROMROOT := proc(poly, r::posint)
local c, exps, i;

    exps := ROOTTOEXPS(r):
    c := poly:

    for i from 1 to nops(exps) do
        c := coeff(c, x[i], exps[i]):
    od:

    return expand(c) mod p;
end:

BUILDSPARSE := proc(coeffs::list, roots::list)
local i, out;

    out := 0:
    for i from 1 to nops(coeffs) do
        out := out + coeffs[i]*ROOTTOMON(roots[i]):
    od:
    return expand(out) mod p:
end:

# ------------------------------------------------------------
# Build moment sequences v_j
# ------------------------------------------------------------
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

printf("Vj's for NUM -> \n"):
vNumer := GETV(ratNumer):
vNumer;

printf("Vj's for DENUM -> \n"):
vDenom := GETV(ratDenum):
vDenom;

# ------------------------------------------------------------
# Hankel matrices for support recovery
# ------------------------------------------------------------
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
local HMat,bVec,L,annP,rootP,R,F,i;

    HMat := paramA:
    bVec := paramB:
    R := table():
    F := table():

    L := LinearSolve(HMat,bVec) mod p:

    annP := [z^tdim,seq(z^(i-1)*L[i],i=1..tdim)]:
    annP := add(annP[i],i=1..nops(annP)):
    annP := Factor(annP) mod p:
    printf("ANNP \n"):
    print(annP);

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

printf("Hankle Matrix for NUM to recover monomials -> \n"):
HNumer := CREATEHANKELMAT(vNumer, termsN):
HNumer;

printf("Hankle Matrix for DENUM to recover monomials -> \n"):
HDenom := CREATEHANKELMAT(vDenom, termsD):
HDenom;

rankHN := Rank(HNumer):
rankHD := Rank(HDenom):
printf("RANK H MAT NUM: %a\nRANK H MAT DENUM: %a\n",rankHN,rankHD):

bVecNumer := GETBVEC(vNumer, termsN):
bVecDenom := GETBVEC(vDenom, termsD):

printf("B vector for NUM -> \n"):
bVecNumer;

printf("B vector for DENOM -> \n"):
bVecDenom;

numRoots,monNumFactor := GETMONFACTOR(HNumer,bVecNumer,termsN):
denomRoots,monDenomFactor := GETMONFACTOR(HDenom,bVecDenom,termsD):

printf("Recovered numerator support factors -> \n"):
monNumFactor;

printf("Recovered denominator support factors -> \n"):
monDenomFactor;

printf("Recovered numerator support count: %a\n", nops(monNumFactor)):
printf("Recovered denominator support count: %a\n", nops(monDenomFactor)):

numRoots := sort(numRoots):
denomRoots := sort(denomRoots):

numCount := nops(numRoots):
denCount := nops(denomRoots):

printf("Recovered numerator roots -> %a\n", numRoots):
printf("Recovered denominator roots -> %a\n", denomRoots):

# ------------------------------------------------------------
# Vandermonde matrices (optional Maple comparison/debug)
# ------------------------------------------------------------
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

printf("VANDERMONDE for NUM -> \n"):
vanMatNumer;

printf("VANDERMONDE for DENOM -> \n"):
vanMatDenom;

GETBVECVAN := proc(paramA::list(integer),n)
local values,bVec,i;

    values := paramA:
    bVec := Vector([seq(values[i],i=1..n)]):
    return bVec:
end proc:

bVecVanNum := GETBVECVAN(vNumer,numCount):
bVecVanDenom := GETBVECVAN(vDenom,denCount):

printf("Vector for numerator Vandermonde solve -> \n"):
bVecVanNum;

printf("Vector for denominator Vandermonde solve -> \n"):
bVecVanDenom;

# ------------------------------------------------------------
# Automatic normalization of the recovered pair
# ------------------------------------------------------------
NORMALIZEPAIR := proc(numCoeff::list, denCoeff::list, denRoots::list, p::prime)
local i, idx, lam, lamInv, newNum, newDen;

    idx := 0:

    # Prefer denominator constant term if present
    for i from 1 to nops(denRoots) do
        if denRoots[i] = 1 then
            idx := i:
            break:
        fi:
    od:

    # Otherwise use first nonzero denominator coefficient
    if idx = 0 then
        for i from 1 to nops(denCoeff) do
            if denCoeff[i] mod p <> 0 then
                idx := i:
                break:
            fi:
        od:
    fi:

    if idx = 0 then
        error "cannot normalize: all denominator coefficients are 0";
    fi:

    lam := denCoeff[idx] mod p:
    if lam = 0 then
        error "cannot normalize: chosen denominator coefficient is 0";
    fi:

    lamInv := 1/lam mod p:

    newNum := [seq(lamInv*numCoeff[i] mod p, i=1..nops(numCoeff))]:
    newDen := [seq(lamInv*denCoeff[i] mod p, i=1..nops(denCoeff))]:

    return newNum, newDen, idx;
end:

NORMALIZEORIGINALPAIR := proc(nPoly, dPoly, denRoots::list, idx::posint, p::prime)
local lam, lamInv;

    lam := COEFFFROMROOT(dPoly, denRoots[idx]) mod p:
    if lam = 0 then
        error "original denominator normalization coefficient is 0";
    fi:

    lamInv := 1/lam mod p:

    return expand(lamInv*nPoly) mod p, expand(lamInv*dPoly) mod p;
end:

# ------------------------------------------------------------
# Coefficient recovery
# VandermondeSolve1 expects (m, v, p, shift)
# ------------------------------------------------------------
LVanNum := VandermondeSolve1([seq(vNumer[i],i=1..numCount)],numRoots,p,1):
LVanDenom := VandermondeSolve1([seq(vDenom[i],i=1..denCount)],denomRoots,p,1):

# Remove the common scalar ambiguity automatically
LVanNum, LVanDenom, normIdx := NORMALIZEPAIR(LVanNum, LVanDenom, denomRoots, p):

recNum := BUILDSPARSE(LVanNum, numRoots):
recDen := BUILDSPARSE(LVanDenom, denomRoots):

# Normalize original pair the same way for exact comparison
nNorm, dNorm := NORMALIZEORIGINALPAIR(n, d, denomRoots, normIdx, p):

printf("RECOVERED NUM AND DENUM (normalized): \n"):
recNum;
recDen;

printf("ORIGINAL NUM AND DENUM with same normalization: \n"):
nNorm;
dNorm;

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

printf("Same rational function? -> %a\n", evalb(expand(recNum*d - recDen*n) mod p = 0)):

printf("Recovered numerator matches normalized original? -> %a\n",
       evalb(expand(recNum - nNorm) mod p = 0)):

printf("Recovered denominator matches normalized original? -> %a\n",
       evalb(expand(recDen - dNorm) mod p = 0)):
