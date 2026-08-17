# This is binding our cpp function into maple. 

addExt := define_external('testAdd64s',
A::ARRAY(0..nn-1,datatype=integer[8]),
B::ARRAY(0..nn-1,datatype=integer[8]),
nn::integer[4],
C::ARRAY(0..nn-1,datatype=integer[8]),
LIB="/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/testAdd.so"):

A := Array(0..2,datatype=integer[8]):
B := Array(0..2,datatype=integer[8]):
C := Array(0..2,datatype=integer[8]): 

A[0]:=1:A[1]:=2:A[2]:=3:
B[0]:=10:B[1]:=20:B[2]:=30:

addExt(A,B,3,C);
[seq(C[i],i=0..2)];

