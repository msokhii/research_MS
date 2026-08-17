lib := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/wrapOBJ2.so":

mVSOLVE := define_external(
    'vSOLVER64C',
    RR::ARRAY(1..nn,datatype=integer[8]),
    LL::ARRAY(1..nn,datatype=integer[8]),
    nn::integer[4],
    aa::ARRAY(1..nn,datatype=integer[8]),
    XX::ARRAY(1..nn,datatype=integer[8]),
    shift::integer[4],
    pp::integer[8],
    LIB=lib
):

VandermondeSolve := proc(m::{Vector,list},v::{Vector,list}, p::prime, shift::integer:=0 )
local t,i,j,M,x,a,q,r,s,temp;
   t := numelems(v);
   if numelems(m) <> t then error "v and m must be the same size"; fi;
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
       for i to t do s := s + v[i]*coeff(q,x,i-1); od;
       a[j] := r*s mod p;
       if shift=0 then next fi;
       r := 1/m[j] mod p;
       r := r &^ shift mod p;
       a[j] := r*a[j] mod p;
   od;
   if type(v,list) then a := convert(a,list); fi;
   return a;
end:

VandermondeSolve1 := proc( v::{Vector,list}, m::{Vector,list}, p::prime, shift::integer:=0 )
local t,i,j,M,x,a,L,R,y;

   t := numelems(v);
   printf("Quadratic C code Vandermonde solver: t=%d  p=%d\n",t,p);
   if numelems(m) <> t then error "v and m must be the same size"; fi;

   R := Array(0..t-1,datatype=integer[8]);
   y := Array(0..t-1,datatype=integer[8]);
   for i from 1 to t do R[i-1] := m[i]; y[i-1] := v[i]; od;

   a := Array(0..t-1,datatype=integer[8]);
   M := Array(0..t,datatype=integer[8]);

   mVSOLVE(R,y,t,a,M,shift,p);
   return [seq( a[i], i=0..t-1 )];

end:

p := 1153;
alpha := numtheory[primroot](p);
e := [11,45,61,72,95,116,201,234,301,310,411,454]:
t := nops(e);
m := [seq( alpha &^ e[i] mod p, i=1..t )]:
c := rand(1000):
f := add( c()*x^e[i], i=1..t ):
f;
v := [seq( Eval(f,x=modp(alpha &^ i,p)) mod p, i=0..t-1 )]:
a := VandermondeSolve( m,v, p ):
a;
add( a[i]*x^e[i], i=1..t ) - f;
a := VandermondeSolve1( m,v, p ):
a;
add( a[i]*x^e[i], i=1..t ) - f;

t := 100:
p := 2^31-1;
alpha := numtheory[primroot](p);
r := rand(1000000):
e := {seq( r(), i=1..t )}:
while nops(e)<t do e := e union {r()} od:
m := [seq( alpha &^ e[i] mod p, i=1..t )]:
c := rand(p):
f := add( c()*x^e[i], i=1..t ):
v := [seq( Eval(f,x=modp(alpha &^ i,p)) mod p, i=0..t-1 )]:
st := time():
a := VandermondeSolve( m,v, p ):
mt := time()-st:
printf("VS1 time=%8.3f\n",mt);
add( a[i]*x^e[i], i=1..t ) - f;
st := time():
a := VandermondeSolve1( m,v, p ):
ft := time()-st:
save a, "aoutput";
printf("VS5 time=%8.3f\n",ft);
add( a[i]*x^e[i], i=1..t ) - f;
