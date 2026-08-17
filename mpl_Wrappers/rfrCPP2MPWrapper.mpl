lib := "/cecm/home/mss59/Desktop/update/pol_ALGO/rfr.so":

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
    LIB=lib
):


# EXAMPLE: 

(* degU := 4: degM := 5:
mLen := degM+1: uLen := degU+1:

M := Array(0..mLen-1, datatype=integer[8], fill=0):
U := Array(0..uLen-1, datatype=integer[8], fill=0):

U[0]:=3: U[1]:=5: U[2]:=7: U[3]:=11: U[4]:=13:
M[5]:=1:

p := 97:
N := 2: DBound := 2:

# OUTPUT BUFFERS: 

nOLEN := N+1:
dOLEN := DBound+1:
nOUT := Array(0..nOLEN-1, datatype=integer[8], fill=0):
dOUT := Array(0..dOLEN-1, datatype=integer[8], fill=0):
degNOUT := 0:
degDOUT := 0:

rc := mRATRECON(mLen,degM,M,uLen,degU,U,N,DBound,p,nOLEN,nOUT,'degNOUT',dOLEN,dOUT,'degDOUT'):
degNOUT;
degDOUT;

print("rc", rc):
print("nOUT full", [seq(nOUT[i], i=0..nOLEN-1)]):
print("dOUT full", [seq(dOUT[i], i=0..dOLEN-1)]):
M;
U; *)

# TESTING PROCEDURE 1: 
amp := 1:
iter2 := 10:
for i from 1 to iter2 do
	d := 100*amp:
	amp := amp*2:
	pp := prevprime(2^63-1):
	printf("DEGREE CHOSEN:%d\nPRIME CHOSEN:%a",d,pp);

	NN := randpoly(x,degree=d,dense) mod pp:
	DD := randpoly(x,degree=d,dense) mod pp:

	while Gcd(NN,DD) mod pp <> 1 do
		DD := randpoly(x,degree=d) mod pp:
	od:

	DD := DD/lcoeff(DD) mod pp:
	MM := mul((x-i),i=1..2*d+1):
	MM := Expand(MM) mod pp:	
	g := Gcdex(DD,MM,x,'SS') mod pp:
	printf("GCD(DD,MM) mod %a = %a\n",pp,g);
	UU := Rem(SS*NN,MM,x) mod pp:
	degU := degree(UU):
	degM := degree(MM):

	M := Array(0..degM, [seq(coeff(MM,x,i),i=0..degM)], datatype=integer[8]):
	U := Array(0..degU, [seq(coeff(UU,x,i),i=0..degU)], datatype=integer[8]):

	# OUTPUT BUFFER INITIALIZATION.
	nOUT := Array(0..degM, datatype=integer[8] ):
	dOUT := Array(0..degM, datatype=integer[8] ):


	rc := mRATRECON(degM+1,degM,M,degU+1,degU,U,d,d,pp,d+1,nOUT,'degNOUT',d+1,dOUT,'degDOUT'):

	# DEGREES: 
	degNOUT;
	degDOUT;

	NNarr := Array(0..degM, [seq(coeff(NN, x, i), i=0..degNOUT)], datatype=integer[8]):
	DDarr := Array(0..degM, [seq(coeff(DD, x, i), i=0..degDOUT)], datatype=integer[8]):
	# MY COMPUTATION: 
	#nOUT; 
	#dOUT;

	# MAPLE: 
	# NN;
	# DD;

	# MY CHECKS SHOULD BE 0.
	nOUT-NNarr;
	dOUT-DDarr;

	#MAPLE CHECK: 
#	mapRat := Ratrecon(UU,MM,x,d,d) mod pp:

	with(Statistics):
	iter := 1:

	TC := Vector(iter):
	TM := Vector(iter):

	for j from 1 to iter do:
        	t := time():
        	repMRC := mRATRECON(degM+1,degM,M,degU+1,degU,U,d,d,pp,d+1,nOUT,'degNOUT',d+1,dOUT,'degDOUT'):
        	TC[j] := time()-t:

#        	t := time():
 #       	repMPLRC := Ratrecon(UU,MM,x,d,d) mod pp:
 #       	TM[j] := time()-t:
	od:

	printf("MY TIMES: \n");
	convert(TC,list);
	printf("MAPLE TIMES: \n");
	convert(TM,list);
	print(d);
	f := fopen("rfrTimingMAPLESERVER.txt",APPEND,TEXT):
	fprintf(f,"DEGREE -> %a\n\n",d);
	fprintf(f,"CPP ROUTINE TIMINGS -> %a\nMAPLE ROUTINE TIMINGS -> %a\n\n",TC,TM);
	fclose(f):
od:

(* with(plots):

dataC := [seq([j, TC[j]], j=1..iter)]:
dataM := [seq([j, TM[j]], j=1..iter)]:

pC := pointplot(dataC, connect=true, symbol=solidcircle):
pM := pointplot(dataM, connect=true, symbol=soliddiamond):

display([pC, pM],
        labels=["trial","seconds"],
               legend=["C wrapper","Maple Ratrecon"]);

Histogram(convert(TC, list), title="C wrapper times"):
Histogram(convert(TM, list), title="Maple Ratrecon times"):*)
