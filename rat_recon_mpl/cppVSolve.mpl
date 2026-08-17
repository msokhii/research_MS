(* CPP version of Vandermonde solver. First run ./runWrap2.sh then compile. *)

libV := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/fastVSolve.so":

mVSOLVE := define_external(
                           'vSOLVER64C',
                           RR::ARRAY(0..nn-1,datatype=integer[8]),
                           LL::ARRAY(0..nn-1,datatype=integer[8]),
                           nn::integer[4],
                           aa::ARRAY(0..nn-1,datatype=integer[8]),
                           XX::ARRAY(0..nn,datatype=integer[8]),
                           shift::integer[4],
                           pp::integer[8],
                           LIB=libV
                           ):

(* Remove type checking from mVSOLVE. *)
mVSOLVE := subsop(1=(
                    RR,
                    LL,
                    nn,
                    aa,
                    XX,
                    shift,
                    pp),
                    op(mVSOLVE)):

cppVS := proc(m,v,p,shift::integer:=0) option inline:
    local t,i,R,y,a,M:

    (* We require V and M to be of the same size. *)
    t := nops(v):
    if nops(m)<>t then
        error "v and m must be the same size":
    fi:

    (* Preparing for input. *)
    R := Array(0..t-1,datatype=integer[8]):
    y := Array(0..t-1,datatype=integer[8]):
    a := Array(0..t-1,datatype=integer[8]):
    M := Array(0..t,datatype=integer[8]):

    for i from 0 to t-1 do
        R[i] := m[i+1]:
        y[i] := v[i+1]:
        a[i] := 0:
    od:
    for i from 0 to t do
        M[i] := 0:
    od:
    mVSOLVE(R,y,t,a,M,shift,p):
    return [seq(a[i], i=0..t-1)]:
end proc: