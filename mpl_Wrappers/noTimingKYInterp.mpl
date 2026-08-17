libV := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/wrapOBJ2.so":
libR := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/wrapOBJ.so":

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

VandermondeSolve1 := proc( m::{Vector,list},v::{Vector,list}, p::prime, shift::integer:=0 )
local t,i,R,y,a,M;

    t := numelems(v);
    if numelems(m) <> t then error "v and m must be the same size"; fi;

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
):

POLYTOARR0 := proc(poly::{polynom,integer}, var::name, deg::integer, p::prime)
local A,i;

    A := Array(0..deg, datatype=integer[8]);
    for i from 0 to deg do
        A[i] := coeff(poly,var,i) mod p;
    od:

    return A;
end:

ARRTOPOLY0 := proc(A::Array, deg::integer, var::name, p::prime)
local i, out;

    if deg < 0 then
        return 0;
    fi;

    out := 0;
    for i from 0 to deg do
        out := out + (A[i] mod p)*var^i;
    od:

    return expand(out) mod p;
end:

Ratrecon1 := proc(Uin::{polynom,integer}, Min::{polynom,integer}, var::name,
                  N::nonnegint, DBound::nonnegint, p::prime)
local Upoly, Mpoly, degU, degM, uLen, mLen,
      UArr, MArr, nOLEN, dOLEN, nOUT, dOUT,
      degNOUT, degDOUT, rc, nn, dd, i;

    Upoly := expand(Uin) mod p:
    Mpoly := expand(Min) mod p:

    degU := degree(Upoly,var):
    degM := degree(Mpoly,var):

    if degU < 0 or degM < 0 then
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

    for i from 0 to nOLEN-1 do nOUT[i] := 0; od:
    for i from 0 to dOLEN-1 do dOUT[i] := 0; od:

    degNOUT := -1:
    degDOUT := -1:

    rc := mRATRECON(
        mLen, degM, MArr,
        uLen, degU, UArr,
        N, DBound, p,
        nOLEN, nOUT, degNOUT,
        dOLEN, dOUT, degDOUT
    ):

    if degNOUT < 0 or degNOUT > N then
        return FAIL;
    fi:

    if degDOUT < 0 or degDOUT > DBound then
        return FAIL;
    fi:

    nn := ARRTOPOLY0(nOUT,degNOUT,var,p):
    dd := ARRTOPOLY0(dOUT,degDOUT,var,p):

    if dd = 0 then
        return FAIL;
    fi:

    return nn/dd;
end:

with(NumberTheory):
p := prevprime(2^31-1):
printf("Prime chosen: %a\n",p):

n := randpoly([x[1],x[2],x[3],x[4]],terms=20) mod p:
d := randpoly([x[1],x[2],x[3],x[4]],terms=20) mod p:
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

(* Assume we also know T >= t where T is the bound on t i.e the actual number of terms.*)
T := max(termsN,termsD):
printf("Value of T: %a\n",T):

(* Assume we know the degrees of both numerator and denominator. *)
N := degree(n):
DD := degree(d):
printf("Degree of NUM: %a\nDegree of DENUM: %a\n",N,DD):

(* To impose uniqueness. *)
test := gcd(n,d):
printf("GCD of NUM and DENUM: %a\n",test):

BLACKBOXF := proc(paramA::list(integer),p)

	local sigmaVal,result,i:
	sigmaVal := paramA:
	result := Eval(f,[seq(x[i]=sigmaVal[i],i=1..nops(indets(f)))]) mod p:
	return result:
end proc:

COMPUTEBETA := proc()

	local betaVal,i,r:
	r := rand(0..p-1):
	betaVal := table():

	for i from 1 to nops(indets(f))-1 do:
    		betaVal[i] := r():
	od:
	betaVal := convert(betaVal,list):
	return betaVal:
end proc:

GETT := proc(j,betaList)

	local sigmaVal,alphaVal,TVal,i,BLACKBOXT,r:
	r := rand(1..p-1):
	TVal := table():	
	sigmaVal := [seq(ithprime(i)^j,i=1..nops(indets(f)))]:
	alphaVal := [seq(r(),i=1..(N+DD+1))]:

    	BLACKBOXT := proc(paramA::integer)

    	local alphaX,result:
    	alphaX := paramA:
    	result := BLACKBOXF([alphaX,seq(betaList[i-1]*alphaX-(betaList[i-1]*sigmaVal[1])+sigmaVal[i],i=2..nops(indets(f)))],p):
    	return result:
    	end proc:

	for i from 1 to (N+DD+1) do:
    		TVal[i] := BLACKBOXT(alphaVal[i]):
	od:

	TVal := convert(TVal,list):
	return (alphaVal,TVal):
end proc:

(* Computing the beta values. *)
betaList := COMPUTEBETA():

INTERSTEP := proc(betaList)

	local alphaVal,TVal,interpVal,ratReconVal,M,i,j:
	ratReconVal := table():
	for j from 1 to 2*(max(termsN,termsD)) do:
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

ratNumer := table():
ratDenum := table():
for i from 1 to nops(ratReconVal) do:
    ratNumer[i] := numer(ratReconVal[i]) mod p:
    ratDenum[i] := denom(ratReconVal[i]) mod p:
od:

printf("Numerators -> \n"):
ratNumer := convert(ratNumer,list):

printf("Denominators -> \n"):
ratDenum := convert(ratDenum,list):

ratNumer:
ratDenum:

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

BUILDSPARSE := proc(coeffs::list, roots::list)
local i, out;

    out := 0:
    for i from 1 to nops(coeffs) do
        out := out + coeffs[i]*ROOTTOMON(roots[i]):
    od:
    return expand(out) mod p:
end:

GETV := proc(paramA::list(polynom))

	local listEval,result,j:
	result := table():
	listEval := paramA:
	for j from 1 to 2*(max(termsN,termsD)) do:
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

with(LinearAlgebra):
CREATEHANKELMAT := proc(paramA::list(integer))

	local val,i,j,H:
	val := paramA:
	H := Matrix(T,T):
	for i from 1 to T do:
    		for j from 1 to T do:
        		H[i,j] := val[i+j-1]:
    		od:
	od:
return H:
end proc:

printf("Hankle Matrix for NUM to recover monomials -> \n"):
HNumer := CREATEHANKELMAT(vNumer):
HNumer;

printf("Hankle Matrix for DENUM to recover monomials -> \n"):
HDenom := CREATEHANKELMAT(vDenom):
HDenom;

(* If T>=t then RANK(H)=t *)
rankHN := Rank(HNumer):
rankHD := Rank(HDenom):
printf("RANK H MAT NUM: %a\nRANK H MAT DENUM: %a\n",rankHN,rankHD):


GETBVEC := proc(paramA::list(integer))

	local values,bVec,i:
	values := paramA:
	bVec := -Vector([seq(values[i],i=T+1..(T+T))]):
	return bVec:
end proc:

bVecNumer := GETBVEC(vNumer):
bVecDenom := GETBVEC(vDenom):

printf("B vector for NUM -> \n"):
bVecNumer;
printf("B vector for DENOM -> \n"):
bVecDenom;

GETMONFACTOR := proc(paramA::Matrix,paramB::Vector)

	local HMat,bVec,L,annP,rootP,R,F,i:

	(* SOLVING HANKEL MATRIX FIRST *)
	HMat := paramA:
	bVec := paramB:
	R := table():
	F := table():
	L := LinearSolve(HMat,bVec) mod p:
	annP := [z^T,seq(z^(i-1)*L[i],i=1..T)]:
	annP := add(annP[i],i=1..nops(annP)):
	annP := Factor(annP) mod p:
	print(annP);
	rootP := Roots(annP) mod p:

	for i from 1 to nops(rootP) do:
    		R[i] := rootP[i][1]:
	od:
	R := convert(R,list):

	for i from 1 to nops(R) do:
    		F[i] := ifactor(R[i]):
	od:
	F := convert(F,list):
	return R,F:
end proc:

numRoots,monNumFactor := GETMONFACTOR(HNumer,bVecNumer):
denomRoots,monDenomFactor := GETMONFACTOR(HDenom,bVecDenom):

monNumFactor;
monDenomFactor;

nops(monNumFactor);
nops(monDenomFactor);

termsN;
termsD;

numRoots := sort(numRoots);
denomRoots := sort(denomRoots);

CREATEVANMAT := proc(paramA::list,n)

	local xVals,i,j,vanMat:
	xVals := paramA:
	vanMat := Matrix(n,n):
	for i from 1 to n do:
    		for j from 1 to n do:
        		vanMat[i,j] := (xVals[j]^(i-1)) mod p:
    		od:
	od:
	
	return vanMat:
end proc:

vanMatNumer := CREATEVANMAT(numRoots,termsN):
vanMatDenom := CREATEVANMAT(denomRoots,termsD):

printf("VANDERMONDE for NUM -> \n"):
vanMatNumer;
printf("VANDERMONDE for DENOM -> \n"):
vanMatDenom;

GETBVECVAN := proc(paramA::list(integer),n)

	local values,bVec,i:
	values := paramA:
	bVec := Vector([seq(values[i],i=1..n)]):
	return bVec:
end proc:

bVecVanNum := GETBVECVAN(vNumer,termsN):
bVecVanDenom := GETBVECVAN(vDenom,termsD):

bVecVanNum;
bVecVanDenom;

LVanNum := VandermondeSolve1([seq(vNumer[i],i=1..termsN)], numRoots, p, 1):
LVanDenom := VandermondeSolve1([seq(vDenom[i],i=1..termsD)], denomRoots, p, 1):

(* FACTOR *)
LI := 1079723992:
LIV := 442395066:

LVanNum := [seq(LIV*LVanNum[i] mod p, i=1..nops(LVanNum))]:
LVanDenom := [seq(LIV*LVanDenom[i] mod p, i=1..nops(LVanDenom))]:

recNum := BUILDSPARSE(LVanNum, numRoots):
recDen := BUILDSPARSE(LVanDenom, denomRoots):

printf("RECOVERED NUM AND DENUM: \n");
recNum;
recDen;

printf("ORIGINAL: \n");
n;
d;

evalb( expand(recNum - n) mod p = 0 );
evalb( expand(recDen - d) mod p = 0 );

(*
scaleDen := [seq(vDenom[j]*(1/trueDen[j] mod p) mod p, j=1..2*T)]:
scaleNum := [seq(vNumer[j]*(1/trueNum[j] mod p) mod p, j=1..2*T)]:

scaleDen;
scaleNum;

convert(scaleDen,set);
convert(scaleNum,set);
*)