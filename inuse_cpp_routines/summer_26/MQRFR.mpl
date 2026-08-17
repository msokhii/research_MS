restart:

p := 101:

HRFTR := proc(M::polynom,UBad::polynom,p::prime,delta::posint)

local R,T,k,LCT,Q,degQ,QMax,NHat,DHat,recF,recG:

R[0],R[1] := M,UBad:
T[0],T[1] := 0,1:
QMax := 0: 
k := 1:

while R[k] <> 0 do
	LCT := (1/lcoeff(T[k],x)) mod p:
	R[k] := (LCT*R[k]) mod p:
	T[k] := (LCT*T[k]) mod p:
	Q := Quo(R[k-1],R[k],x) mod p:
	degQ := degree(Q):
        print(degQ):
	if degQ > QMax then
		QMax := degQ:
		NHat := R[k]:
		DHat := T[k]:
        elif QMax = 2 and degQ >= QMax then
		print(n,d,R[k],T[k]):
	fi:
	R[k+1] := Expand(R[k-1]-(Q*R[k])) mod p:
	T[k+1] := Expand(T[k-1]-(Q*T[k])) mod p:
	k := k+1:
od: 
if QMax < 2 then
	return FAIL:
fi:

Lambda := Gcd(NHat,DHat) mod p:
recF := Quo(NHat,Lambda,x) mod p:
recG := Quo(DHat,Lambda,x) mod p:
print(Lambda);
return (recF/recG):
end proc:

MQRFR := proc(M::polynom,U::polynom,p::prime)

local R,T,k,LCT,Q,degQ,QMax,n,d:

R[0],R[1] := M,U: 
T[0],T[1] := 0,1:
k := 1:
QMax := 0:

while R[k] <> 0 do
	LCT := (1/lcoeff(T[k],x)) mod p:
	R[k] := (LCT*R[k]) mod p:
	T[k] := (LCT*T[k]) mod p:
	Q := Quo(R[k-1],R[k],x) mod p:
        degQ := degree(Q):
	if degQ > QMax then
		QMax := degQ:
		n := R[k]:
		d := T[k]:
	elif QMax = 2 and degQ >= QMax then
		print(n/d,R[k],T[k]):
	fi:
	R[k+1] := Expand(R[k-1]-(Q*R[k])) mod p:
	T[k+1] := Expand(T[k-1]-(Q*T[k])) mod p:
	k := k+1:
od:

if QMax < 2 then
	return FAIL:
fi:
return (n/d):
end proc:

n := (x^2+3*x+1) mod p:
d := (x+2) mod p:
degN := degree(n):
degD := degree(d): 

(* For MQRFR we need degN+degD+2 points. *)
pts := degN+degD+2:

getUandM := proc(n::polynom,d::polynom,p::prime,pts::posint)

local U,M,i,y,evalD,evalN,g:

y := table():
for i from 1 to pts do
	evalD := Eval(d,x=i) mod p:
	if evalD = 0 then
		return FAIL:
	fi:
	evalN := Eval(n,x=i) mod p:
	g := Gcdex(evalD,p,x,'s','t') mod p:
        if g<> 1 then
		return FAIL:
	fi:
	y[i] := Expand(evalN*s) mod p: 
od:

y := convert(y,list):
U := Interp([seq(i,i=1..pts)],y,x) mod p:
M := Expand(mul(x-i,i=1..pts)) mod p: 

return (M,U):
end proc:

getUandM2 := proc(n::polynom,d::polynom,p::prime,pts2::posint)

local U,M,i,y,evalD,evalN,g:

y := table():
for i from 1 to pts2 do
	evalD := Eval(d,x=i) mod p:
	if evalD = 0 then
		return FAIL:
	fi:
	evalN := Eval(n,x=i) mod p:
	g := Gcdex(evalD,p,x,'s','t') mod p:
	if g <> 1 then
		return FAIL:
	fi:
	y[i] := Expand(evalN*s) mod p:
od:

y := convert(y,list):
yBad := subsop(2=(y[2]+10) mod p,y):
print(y):
print(yBad):
U := Interp([seq(i,i=1..pts2)],yBad,x) mod p:
M := Expand(mul(x-i,i=1..pts2)) mod p:

return (M,U):
end proc:

M,U := getUandM(n,d,p,pts):
M;
U;
M2,U2 := getUandM2(n,d,p,degN+degD+2+2):
M2;
U2;

recover := MQRFR(M,U,p):
recover;

check := ((n/d)-recover) mod p:
check;

recover2 := HRFTR(M2,U2,p,5):
recover2;
