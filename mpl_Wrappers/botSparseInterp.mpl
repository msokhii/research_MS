restart:
with(LinearAlgebra): 
M1 := x[1]^3*x[2]^4:
M2 := x[1]*x[2]^3*x[3]:
M3 := x[1]^6*x[3]^2: 
a1 := 101:
a2 := 103: 
a3 := 105:
(* 
f belongs in Z[x[1],x[2],x[3]] and this is precisely what we want to recover. We want to recover the 
coefficients (a[1],a[2],a[3]) and the monomials efficiently. 
*)

f := a1*M1+a2*M2+a3*M3:
print(f):

(*
We will make use of Ben-Or Tiwari sparse interpolation technique. This algorithm requires 2T points
from Z, where T>=t and t is number of terms in our polynomial. In our f, t=3. So lets, pick T=4. Clearly, 
T>=t so this method is guaranteed to work.

STEP 1: 

We will compute vj = f(p1^j,p2^j,p3^j) for distinct primes p1,p2,p3 and 0<=j<=2T-1.
*)
computeVJ := proc(f::polynom,T::posint)

local j,res:
res := table():
for j from 0 to 2*T-1 do
    res[j] := eval(f,{x[1]=2^j,x[2]=3^j,x[3]=5^j}):
od:
res := convert(res,list):
return res:
end proc:
vVal := computeVJ(f,4):
vVal;
nops(vVal);

(* 
Recall, the manipulation we did for the lambda polynomial. This gives rise to a T x T Hankel matrix. We are 
solving the linear system  for lambda H*lambda=v.
*)  
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
H := createHMatrix(vVal,4):
H;
      
(*
Theorem:
If T>=t, then rank(H)=t. We can use this theorem to check that we indeed have the correct value for T.
*)

compRank := Rank(H):
compRank;
                              
H := H[1..3,1..3]:
H;
     
createVHVector := proc(vVal::list(integer),T::posint,compRank::posint)

local v,i:
v := Vector[column](compRank):
for i from 1 to compRank do
    v[i] := -vVal[i+T-1]:
od:
return v:
end proc:
v := createVHVector(vVal,4,3):
v;
        
(*
These values are L0,L1 and L2 respectively.
*)

L := LinearSolve(H,v):
L;
     
(*
We can now construct the Lambda polynomial i.e. sum from j=0 to t.
*) 

Lambda := z^3+L[3]*z^2+L[2]*z+L[1]:
Lambda;
          
(* 
We now factor the Lambda polynomial and get the roots.
*)

factor(Lambda);
         
R := roots(Lambda):
R;
           
m1,m2,m3 := seq(r[1],r in R):
ifactor(m1);
ifactor(m2);
ifactor(m3);
                
(*
We have succesfully recovered the monomials.
*)

mon1 := x[1]^6*x[3]^2:
mon2 := x[1]^3*x[2]^4:
mon3 := x[1]*x[2]^3*x[3]:
print(mon1):
print(mon2):
print(mon3):
        
(* 
Now, we want to recover the unknown coefficients ai's. Recall, this gives rise to a transposed Vandermonde
matrix.
*) 

V := Matrix([[1,1,1],[m1,m2,m3],[m1^2,m2^2,m3^2]]):
V;
          
createSmallV := proc(vVal::list(integer),compRank::posint)

local v,i:
v := Vector[column](compRank):
for i from 1 to compRank do
    v[i] := vVal[i]:
od:
return v:
end proc:
smallV := createSmallV(vVal,compRank):
smallV;
     
(*
We now solve the linear system for a i.e. Va=v.
*)

A := LinearSolve(V,smallV):
A;
     
recoveredF := A[1]*mon1+A[2]*mon2+A[3]*mon3:
print(recoveredF):
       
print(f):
