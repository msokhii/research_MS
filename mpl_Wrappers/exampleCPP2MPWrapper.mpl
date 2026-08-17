lib := "/local-scratch/localhome/mss59/Desktop/research_Works/development/libpoly.so":
pADDIP64 := define_external(
  'pADDIP64_C',
  a_len::integer[4],
  A::ARRAY(0..a_len-1, datatype=integer[8]),
  degA::integer[4],
  b_len::integer[4],
  B::ARRAY(0..b_len-1, datatype=integer[8]),
  degB::integer[4],
  p::integer[8],
  RETURN::integer[4],
  LIB=lib
):

# Example polys: a(x)=3+2x+5x^2+4x^3, b(x)=1+x^3
degA := 3: degB := 3:
maxDeg := max(degA, degB):
a_len := maxDeg+1:
b_len := degB+1:

A := Array(0..a_len-1, datatype=integer[8], fill=0):
B := Array(0..b_len-1, datatype=integer[8], fill=0):

A[0]:=3: A[1]:=2: A[2]:=5: A[3]:=4:
B[0]:=1: B[3]:=1:

p := 17:

newDeg := pADDIP64(a_len, A, degA, b_len, B, degB, p):
print(newDeg);

# Result coefficients are now in A[0..newDeg]
seq(A[i], i=0..newDeg);

