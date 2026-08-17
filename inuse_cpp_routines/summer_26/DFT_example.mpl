restart: 

p := 101:
f := Randpoly(3,x) mod p:
g := Randpoly(1,x) mod p: 
checkGCD := Gcd(f,g) mod p:

f;
g;
checkGCD;

LCG := lcoeff(g,x) mod p:
g := Expand((1/LCG)*g) mod p:
checkGCD := Gcd(f,g) mod p:

g;
checkGCD;

h := f/g:
h;

degF := degree(f,x):
degG := degree(g,x):

degF;
degG;

E := 2: 
pts := degF+degG+(2*E)+1:

E;
pts;

getUandM := proc(pts::posint,p::prime)

local numEval,denEval,yVal,i,k,M,U: 

yVal := table():
for i from 1 to pts do
	denEval := Eval(g,x=i) mod p:
	if denEval = 0 then
		return FAIL:
	fi:
	numEval := Eval(f,x=i) mod p:
	yVal[i] := numEval/denEval mod p:
od:

yVal := convert(yVal,list):
M := Expand(mul(x-i,i=1..pts)) mod p:
U := Interp([seq(i,i=1..pts)],yVal,x) mod p:

return (M,U,yVal):
end proc:

M,U,yVal := getUandM(pts,p):

M;
U;
yVal;

yValC := [91, 69, 1, 65, 85, 46, 55, 57, 22]:
yValC; 

UC := Interp([seq(i,i=1..pts)],yValC,x) mod p:
UC;

DFTRFR := proc(M,U,degF::integer,degG::integer,E::integer)

local degM,checkGCD,degGCD,R,T,Q,d,p2,q,k,LCQ: 

degM := degree(M,x):
if degM<degF+degG+(2*E) then
	return FAIL:
fi:

checkGCD := Gcd(M,U) mod p:
degGCD := degree(checkGCD,x):

if degGCD>degF+E then
	return FAIL:
fi: 

R[0] := M:
R[1] := U: 
T[0] := 0:
T[1] := 1:
k := 1: 

while R[k]<>0 and degree(T[k],x)<=degG+E do
	Q[k] := Quo(R[k-1],R[k],x) mod p:
	R[k-1] := R[k]:
	R[k] := Expand(R[k-1]-Q[k]*R[k]) mod p:
        T[k-1] := T[k]:
	T[k] := Expand(T[k-1]-Q[k]*T[k]) mod p:
        k := k+1:
od: 

d := Gcd(R[k-1],T[k-1]) mod p:
p2 := Expand(R[k-1]/d) mod p: 
q := Expand(T[k-1]/d) mod p: 

LCQ := lcoeff(q,x) mod p:
p2 := Expand((1/LCQ)*p) mod p:
q := Expand((1/LCQ)*q) mod p:

return (p2,q):
end proc:


DFTRFR(M,UC,degF,degG,E);
