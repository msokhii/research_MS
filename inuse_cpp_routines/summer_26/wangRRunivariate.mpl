restart:

p := 2^31-1: 
r := rand(p):
degN := 80:
degD := 80: 
n := randpoly(x,coeffs=r,degree=degN): 
d := randpoly(x,coeffs=r,degree=degD):
g := Gcd(n,d) mod p:

while g <> 1 do
	n := randpoly(x,coeffs=r,degree=degN):
	d := randpoly(x,coeffs=r,degree=degD):
        g := Gcd(n,d) mod p:
od:

normD := proc(n,d)
local LC,g,inv,N2,D2: 

N2 := n:
D2 := d:
LC := lcoeff(D2): 
if LC<>1 then
	g := Gcdex(LC,p,x,'s','t') mod p:
	if g<>1 then
		return 'FAIL':
	fi:
	inv := s:
	N2 := Expand(N2*inv) mod p:
	D2 := Expand(D2*inv) mod p:
fi:

return (N2,D2):
end proc: 

n,d := normD(n,d):
n;
d;

(* ATP we have LC(d)=1 and GCD(n,d)=1 to impose uniqueness. *) 
(* Since, we know the degrees of n and d beforehand we know exactly the number of points needed for recovery. *) 

m := degN+degD+1: 
alphaVal := Array(1..m):
for i from 1 to m do
	alphaVal[i] := i: 
od:

evalND := proc(n,d)
local yVal,evalD,evalN,g,inv,i,N,D: 

N := n:
D := d:
yVal := Array(1..m):
for i from 1 to m do
	evalD := Eval(D,x=alphaVal[i]) mod p:
	if evalD=0 then
		return 'FAIL':
	fi:
	evalN := Eval(N,x=alphaVal[i]) mod p:
	g := Gcdex(evalD,p,x,'s','t') mod p:
	if g<>1 then
		return 'FAIL':
	fi:
	inv := s:
	yVal[i] := Expand(evalN*inv) mod p:
od:

return yVal:
end proc:

yVal := evalND(n,d): 
yVal;

intTime := time():
to 10^5 do:
U := Interp(alphaVal,yVal,x) mod p:
od:
intStop := time()-intTime:
printf("Avg call time for interp: %a\n",intStop/10^5):
M := Expand(mul(x-alphaVal[i],i=1..m)) mod p:

U;
M; 

degM := degree(M):
g := Gcd(d,M) mod p:
g;
degM; 

(* Now, GCD(D,M)=1 and deg(M)=deg(N)+deg(D)+1 -> We can use Wang's RFR. *) 

EEA := proc(U,M) 
local r,t,u,k,res,remSeq,quoSeq,tSeq: 

r[0] := M:
r[1] := U: 
t[0] := 0:
t[1] := 1: 
remSeq := table():
quoSeq := table():
tSeq := table():
k := 1: 

while r[k]<>0 do 
	u := 1/lcoeff(t[k]) mod p:
        r[k] := u*r[k] mod p:
	t[k] := u*t[k] mod p:
	remSeq[k] := degree(r[k]):
	tSeq[k] := degree(t[k]):
	q := Quo(r[k-1],r[k],x) mod p:
        quoSeq[k] := degree(q[k]):
	print(Step(k),r[k]/t[k],Quo(q)):
	r[k+1] := Expand(r[k-1]-q*r[k]) mod p:
	print(r[k+1]):
	t[k+1] := Expand(t[k-1]-q*t[k]) mod p:
	print(t[k+1]):
	k++:
od: 

remSeq := convert(remSeq,list):
quoSeq := convert(quoSeq,list):
tSeq := convert(tSeq,list): 
print(remSeq):
print(quoSeq):
print(tSeq):
res := r[k-1]:
return res:
end proc:

compRes := EEA(U,M);
f := n/d:
f;
f-compRes;
