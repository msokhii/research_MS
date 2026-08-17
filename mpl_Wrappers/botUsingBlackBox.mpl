restart:
with(LinearAlgebra):
p := prevprime(2^32-1):
print(p):
f := randpoly([x[1],x[2],x[3]],degree=10):
print(f):
nops(f);

(* Is this sparse? Recall, for a polynomial f to be sparse we require 
   t <= sqroot(product from 1 to n of (di+1)). *)


checkSparse := proc(f::polynom)

local terms,t,i:
t := 1:
terms := nops(f):
for i from 1 to nops(indets(f)) do
    t := t*(degree(f,x[i])+1):
od:
t := evalf(sqrt(t)):
if terms<=t then
    return 1:
else
    return 0:
fi:
end proc:
(* Clearly, f is sparse and we can proceed. *)

check1 := checkSparse(f):
check1;
          
BF := proc(val::list(integer),p::prime)

local res,i:
res := Eval(f,[seq(x[i]=val[i],i=1..nops(indets(f)))]) mod p:
return res:
end proc:
(* Assuming, we know t. Let T>=t i.e. T=6. *)

computeVvals := proc(T::posint)

local vVal,i:
vVal := table():
for i from 0 to (2*T-1) do
    vVal[i] := BF([2^i,3^i,5^i],p): 
od:
vVal := convert(vVal,list):
return vVal:
end proc:
vVal := computeVvals(6):
vVal;
nops(vVal);

(* We will first find the monomials. *)
createHMatrix := proc(vVal::list(integer),T::posint)

local H,i,j: 
H := Matrix(T,T):
for i from 1 to T do:
    for j from i to T do:
        (* As, H is square. We can have some speedup. *)
        H[i,j] := vVal[i+j-1]:
        H[j,i] := H[i,j]:
    od:
od:
return H:
end proc:
H := createHMatrix(vVal,6):
H;
    
(* We can invoke the following theorem: If T>=t then rank(H)=t. *)

compRank := Rank(H):
compRank;
                              
createVHVector := proc(vVal::list(integer),T::posint,compRank::posint)

local v,i:
v := Vector[column](compRank):
for i from 1 to compRank do
    v[i] := -vVal[i+T]:
od:
return v:
end proc:
v := createVHVector(vVal,6,6):
v;
  
(* Now we can solve the linear system: H*Lambda=v. *)

L := LinearSolve(H,v) mod p:
L;
      
(*
We can now construct the Lambda polynomial i.e. sum from j=0 to t.
*) 

Lambda := z^6+L[6]*z^5+L[5]*z^4+L[4]*z^3+L[3]*z^2+L[2]*z+L[1]:
Lambda;

(* We can now factor the lambda polynomial to get the roots. *)

Factor(Lambda) mod p;
R := Roots(Lambda) mod p:
R;

m1,m2,m3,m4,m5,m6 := seq(r[1],r in R):
ifactor(m1);
ifactor(m2);
ifactor(m3);
ifactor(m4);
ifactor(m5);
ifactor(m6);
           
mon1 := x[1]*x[2]:
mon2 := x[1]^4*x[3]:
mon3 := x[1]^4*x[2]*x[3]^5:
mon4 := x[1]^2*x[2]^5:
mon5 := x[1]^2*x[2]^4*x[3]^3:
mon6 := x[1]^5*x[2]^2*x[3]^3:
print(mon1):
print(mon2):
print(mon3):
print(mon4):
print(mon5):
print(mon6):
             
(* We have succesfully recovered the monomials. Now we want to recover the coefficients. But recall, 
to recover the coefficients, this gives rise to a transposed vandermonde matrix. *)
V := Matrix([[1,1,1,1,1,1],[m1,m2,m3,m4,m5,m6],[m1^2,m2^2,m3^2,m4^2,m5^2,m6^2],
             [m1^3,m2^3,m3^3,m4^3,m5^3,m6^3],[m1^4,m2^4,m3^4,m4^4,m5^4,m6^4],
             [m1^5,m2^5,m3^5,m4^5,m5^5,m6^5]]) mod p:
V;
          
createSmallV := proc(vVal::list(integer),compRank::posint)

local v,i:
v := Vector[column](compRank):
for i from 1 to compRank do
    v[i] := vVal[i] mod p:
od:
return v:
end proc:
smallV := createSmallV(vVal,compRank):
smallV;
      
(*
We now solve the linear system for a i.e. Va=v.
*)

A := LinearSolve(V,smallV) mod p:
A;
   
f;

