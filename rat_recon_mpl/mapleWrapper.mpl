libObj := "/cecm/home/mss59/Desktop/research_ServerMaple/inuse_cpp_routines/cpp_routines/cppObj_old.so":
(* libObj := "/Users/msokhi/Desktop/researchFiles/newDir/routinesCPP/cppObj.so": *)

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
                            LIB=libObj
                            ): 

(* Remove type checking from mRATRECON. *)

mRATRECON := subsop(1=(
                       mLen,
                       degM,
                       M,
                       uLen,
                       degU,
                       U,
                       N,
                       DBound,
                       p,
                       nOLEN,
                       nOUT,
                       degNOUT,
                       dOLEN,
                       dOUT,
                       degDOUT),
                       op(mRATRECON)):

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
                                LIB=libObj
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

mVSOLVE := define_external('cppVSolve',
  mLen::integer[4],
  mIn::ARRAY(0..mLen-1,datatype=integer[8]),
  yLen::integer[4],
  yIn::ARRAY(0..yLen-1,datatype=integer[8]),
  shiftInt::integer[4],
  pp::integer[8],
  outLen::integer[4],
  aOut::ARRAY(1..outLen,datatype=integer[8]),
  RETURN::integer[4],
  LIB=libObj):

(* Remove typechecking from mVSOLVE *)

mVSOLVE := subsop(1=(
                    mLen,
                    mIn,
                    yLen,
                    yIn,
                    shiftInt,
                    pp,
                    outLen,
                    aOut),op(mVSOLVE)): 

mBM := define_external('cppBM',
  aLen::integer[4],
  aIn::ARRAY(0..aLen-1, datatype=integer[8]),
  p::integer[8],
  outLen::integer[4],
  lOut::ARRAY(0..outLen-1, datatype=integer[8]),
  degOut::REF(integer[4]),
  RETURN::integer[4],
  LIB=libObj):

mBM := subsop(1=(
                aLen,
                aIn,
                p,
                outLen,
                lOut,
                degOut),op(mBM)): 

(* Maple wrapper for Berlekamp Massey *)

cppBMM := proc(a::{Vector,list},p::prime) option inline:
local N,n,i,aArr,L,deg,rc:

   N := numelems(a):
   if N <= 0 then 
        error "Input sequence must be nonempty.":
   fi:

   aArr := Array(0..N-1,datatype=integer[8]):
   for i from 1 to N do 
        aArr[i-1] := a[i]: 
   od:

   n := iquo(N,2):
   L := Array(0..n,datatype=integer[8]):
   deg := -1:
   
   rc := mBM(N,aArr,p,n+1,L,deg):
   if rc <> 0 then 
        error "BM cpp version returned error code %1", rc: 
   fi:
   
   deg := -1:
   for i from n by -1 to 0 do
       if L[i] <> 0 then 
            deg := i:
            break
       fi:
   od:
  
   if deg < 0 then
        return []: 
   fi:
   return [seq(L[i],i=0..deg)]:
end proc:

(* Maple wrapper for Vandermonde Solve. *)

cppVS := proc(v::{Vector,list},m::{Vector,list},p::prime,shift::integer:=0) option inline:
local t,i,a,R,y,rc:

   t := numelems(v):
   if numelems(m) <> t then 
        error "v and m must be the same size.": 
   fi:

   R := Array(0..t-1,datatype=integer[8]):
   y := Array(0..t-1,datatype=integer[8]):
   for i from 1 to t do 
        R[i-1] := m[i]:
        y[i-1] := v[i]: 
   od:

   a := Array(0..t-1,datatype=integer[8]):

   rc := mVSOLVE(t,R,t,y,shift,p,t,a):
   if rc <> 0 then 
        error "VSolve cpp version returned error code %1", rc:
   fi:

   return [seq(a[i],i=0..t-1)];
end proc:

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

(* Converts a polynomial in Maple rep. to an array of coeffs. *)
convertPY2ARR := proc(poly,var,deg,p) option inline:
    local A,i;
    
    A := Array(0..deg,datatype=integer[8]):
    for i from 0 to deg do
        A[i] := coeff(poly,var,i):
    od:
    return A:
end proc:

(* Converts an array of coeffs. to a polynomial in Maple rep. 
   Array coeffs. are from low deg to high. *)
convertARR2PY := proc(A,deg,var) option inline:
    local i:
    
    add(A[i]*var^i,i=0..deg):
end proc:

(* Check if deg(A[i])=0 *)
checkZeroPY := proc(A,len) option inline:
    local i:

    for i from len-1 by -1 to 0 do
        if A[i]<>0 then
            return i:
        fi:
    od:
    print("ZERO POLYNOMIAL."):
    return -1:
end proc:

lastOP := 0:

(* Maple wrapper for mRATRECON. We call this maple function. *)
cppRR := proc(Uin,
              Min, 
              var,
              N, 
              DBound,
              p) option inline:
    
    global lastOP:
    local Upoly,Mpoly,degU,degM,uLen,mLen,
          UArr,MArr,nOLEN,dOLEN,nOUT,dOUT,
          degNOUT,degDOUT,cppRet,nn,dd,i:

        Upoly := Uin: #U
        Mpoly := Min: #M
        degU := degree(Upoly,var):
        degM := degree(Mpoly,var):

        if degU<0 or degM<0 then
            lastOP := -999: #U or M is a 0 polynomial.
            error "U or M is a 0 polynomial.":
        fi:

        uLen := degU+1:
        mLen := degM+1:

        (* Prepare inputs for cpp routine. *)
        UArr := convertPY2ARR(Upoly,var,degU,p):
        MArr := convertPY2ARR(Mpoly,var,degM,p):

        (* This is what the cpp routine will output. *)
        nOLEN := N+1:
        dOLEN := DBound+1:
        nOUT := Array(0..nOLEN-1,datatype=integer[8]):
        dOUT := Array(0..dOLEN-1,datatype=integer[8]):

        (* Initial State. Assuming N and D output will not be the 0 polynomial. *)
        degNOUT := -1:
        degDOUT := -1:
         
        cppRet := mRATRECON(
                            mLen,degM,MArr,
                            uLen,degU,UArr,
                            N,DBound,p,
                            nOLEN,nOUT,degNOUT,
                            dOLEN,dOUT,degDOUT
        ):
        lastOP := cppRet:
        
        (* 0 flag means success in reconstruction. *)
        if cppRet <> 0 then
            return FAIL:
        fi:
        if degNOUT<0 or degNOUT>N then
            degNOUT := checkZeroPY(nOUT,nOLEN);
        fi:
        if degDOUT<0 or degDOUT>DBound then
            degDOUT := checkZeroPY(dOUT,dOLEN);
        fi:
        if degNOUT<0 or degNOUT>N then
            error "degNOUT<0 or degNOUT>N.":
        fi:
        if degDOUT<0 or degDOUT>DBound then
            error "degDOUT<0 or degDOUT>DBound.":
        fi:
        
        (* Prepare output for maple rep. *)
        nn := convertARR2PY(nOUT,degNOUT,var):
        dd := convertARR2PY(dOUT,degDOUT,var):

        if dd=0 then
            error "DD=0.":
        fi:

        (* Return reconstruction. *)
        return (nn/dd):
end proc:

(* 
This version takes in a list instead of an array. 
*)
cppNewtonInterp := proc(xVals,yVals,var,p) option inline:
    global lastOP:
    local n,xArr,yArr,outLen,yOut,degOut,cppRet,poly,i:

    n := nops(xVals):

    (* We want f(x_i)=y_i so these should have the same size. *)
    if n<>nops(yVals) or n=0 then
        lastOP := -100: #Flag for not the same size.
        error "X and Y are not the same size.":
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
    fi:
    degOut := checkZeroPY(yOut,outLen):
    if degOut<0 then
        error "0 polynomial.":
    fi:
    poly := add(yOut[i]*var^i,i=0..degOut) mod p:
    return poly:
end proc:
