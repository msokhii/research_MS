restart:
with(NumberTheory):
p := prevprime(2^63-1): 

g := randpoly([x[1],x[2]],dense,degree=5,terms=10) mod p:
h := randpoly([x[1],x[2]],dense,degree=3,terms=5) mod p:

while gcd(g,h) <> 1 do
	h := randpoly([x[1],x[2]],dense,degree=3,terms=5):
od:

printf("G: %a\nH: %a\n",g,h):

f := g/h:

printf("F: %a\n",f):

# Assume, we know the degrees of g(numerator) and h(denominator). 

n := degree(g):
dd := degree(h): 

printf("DEGREE of G: %a\nDEGREE OF H: %a\n",n,dd):

# Only, purpose of this is to bound j.

termsG := nops(g):
termsH := nops(h):

printf("TERMS in G: %a\nTERMS in H: %a\n",termsG,termsH):

printf("PRIME CHOSEN: %a\n",p):

# Assume, we get f from a modular black box B: F^n -> F. 

BB := proc(paramA::list(integer),p)

local sigmaVal,result,i: 

sigmaVal := paramA:
result := Eval(f,[seq(x[i]=sigmaVal[i],i=1..nops(indets(f)))]) mod p:
return result:
end proc:

# KY, reduce the problem from interpolating f directly to interpolating 
# u*g and u*h from some unit u in F. For this, they define the following mapping: x -> x_1 and
# x_i -> B_1*x_i - B_i(sigma_1-sigma_i).

# We need B_2,B_3, ... , B_n values at random. 

getBeta := proc()

local val,i,r: 

r := rand(0..p-1):
val := table():
for i from 1 to nops(indets(f))-1 do
	val[i] := r():
od:
val := convert(val,list):
return val:
end proc:

TVal := proc(j,betaVal)

local sigmaVal,alphaVal,tVal,i,BBT,r: 

r := rand(1..p-1):
tVal := table():
sigmaVal := [seq(ithprime(i)^j,i=1..nops(indets(f)))]:
alphaVal := [seq(r(),i=1..(n+dd+1))]:

printf("SIGMA VALUES: %a\nALPHA VALUES: %a\n",sigmaVal,alphaVal): 

	BBT := proc(paramA::integer)
	
	local temp,res:

	temp := paramA:
	res := BB([temp,seq(betaVal[i-1]*temp-((betaVal[i-1]*sigmaVal[1])+sigmaVal[i]),i=2..nops(indets(f)))],p):
	return res:
	end proc:

for i from 1 to (n+dd+1) do
	tVal[i] := BBT(alphaVal[i]):
od:

tVal := convert(tVal,list):
return(alphaVal,tVal):
end proc: 

betaVal := getBeta(): 

printf("BETA VALUES: %a\n",betaVal):

STEP2 := proc(betaVal)

local j,alphaVal,tVal,interpVal,M,ratCon:
 
ratCon := table():
for j from 1 to 2*(max(termsG,termsH)) do: 
	alphaVal,tVal := TVal(j,betaVal):
	interpVal := Interp(alphaVal,tVal,z) mod p:
	M := [seq(z-alphaVal[i],i=1..nops(alphaVal))]:
	M := Expand(convert(M,`*`)) mod p:
	ratCon[j] := RatRecon(interpVal,M,z,n,dd) mod p:
od:
ratCon := convert(ratCon,list):
return ratCon:
end proc:

ratVal := STEP2(betaVal): 
ratVal;

ratN := table():
ratD := table(): 

for i from 1 to nops(ratVal) do
	ratN[i] := numer(ratVal[i]): 
	ratD[i] := denom(ratVal[j]): 
od: 

ratN := convert(ratN,list):
ratD := convert(ratD,list): 

ratN;
ratD;
