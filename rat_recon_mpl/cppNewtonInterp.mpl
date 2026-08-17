(* CPP version of Maple's Interp. First run ./runNewtonInterp.sh then compile. *)

libNewton := "/local-scratch/localhome/mss59/Desktop/research_Works/development/pol_ALGO/newton.so":

mNEWTONINTERP := define_external(
                                'cppInterp',
                                xLen::integer[4],
                                xIn::ARRAY(0..xLen-1,datatype=integer[8]),
                                yLen::integer[4],
                                yIn::ARRAY(0..yLen-1,datatype=integer[8]),
                                p::integer[8],
                                outLen::integer[4],
                                yOut::ARRAY(0..outLen-1,datatype=integer[8]),
                                degOut::REF(integer[4]),
                                RETURN::integer[4],
                                LIB=libNewton
                                ):

(* Remove type checking from mNEWTONINTERP. *)

mNEWTONINTERP := subsop(1=(
                           xLen,
                           xIn,
                           yLen,
                           yIn,
                           p,
                           outLen,
                           yOut,
                           degOut),
                           op(mNEWTONINTERP)):

(* Check if deg(A[i])=0 *)
checkZeroPY := proc(A,len) option inline:
    local i:

    for i from len-1 by -1 to 0 do
        if A[i]<>0 then
            return i:
        fi:
    od:
    return -1:
end proc:

lastOP := 0:

cppNewtonInterp := proc(xVals,yVals,var,p) option inline:
    global lastOP;
    local n,xArr,yArr,outLen,yOut,degOut,cppRet,poly,i;

    n := nops(xVals):

    (* We want f(x_i)=y_i so these should have the same size. *)
    if n<>nops(yVals) or n=0 then
        lastOP := -100: #Flag for not the same size.
        print("X and Y are not the same size : ",lastOP):
        return FAIL:
    fi:

    (* Preparing for input. *)
    xArr := Array(0..n-1,datatype=integer[8]):
    yArr := Array(0..n-1,datatype=integer[8]):

    for i from 1 to n do
        xArr[i-1] := xVals[i]:
        yArr[i-1] := yVals[i]:
    od:

    outLen := n:
    yOut := Array(0..outLen-1,datatype=integer[8]):
    degOut := 0:

    cppRet := mNEWTONINTERP(
                            n,xArr,
                            n,yArr,
                            p,
                            outLen,yOut,degOut
    ):

    lastOP := cppRet:
    if cppRet <> 0 then
        return FAIL:
    else
        print("INTERPOLATION SUCCESFUL!"):
    fi:
    degOut := checkZeroPY(yOut,outLen):
    if degOut<0 then
        print("ZERO POLYNOMIAL : ",degOut):
        return FAIL:
    fi:
    poly := add(yOut[i]*var^i,i=0..degOut) mod p:
    return poly:
end proc:
