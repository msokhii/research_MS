 libObj := "/cecm/home/mss59/Desktop/resMaple_MS/inuse_cpp_routines/cpp_routines/cppObj.so":
#libObj := "/home/msokhi/Desktop/res_MS/inuse_cpp_routines/cpp_routines/cppObj.so":

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

mINTERPDFTRFR := define_external(
                          'cppInterpDFTRFR',
                          nPts::integer[4],
                          alpha::ARRAY(0..nPts-1, datatype=integer[8]),
                          Yin::ARRAY(0..nPts-1, datatype=integer[8]),
                          N::integer[4],
                          DBound::integer[4],
                          E::integer[4],
                          p::integer[8],
                          nOLEN::integer[4],
                          nOUT::ARRAY(0..nOLEN-1, datatype=integer[8]),
                          dOLEN::integer[4],
                          dOUT::ARRAY(0..dOLEN-1, datatype=integer[8]),
                          RETURN::integer[4],
                          LIB=libObj
                          ):

mINTERPDFTRFR := subsop(1=(
                    nPts,alpha,Yin,
                    N,DBound,E,p,
                    nOLEN,nOUT,
                    dOLEN,dOUT),
                    op(mINTERPDFTRFR)):


cppInterpDFTRFR := proc(xVals, yVals, var, N, DBound, E, p)
    global lastOP:
    local n,xArr,yArr,nOLEN,dOLEN,nOUT,dOUT,cppRet,degN,degD,i,nn,dd:

        n := nops(xVals):
        if n <> nops(yVals) or n=0 then
            lastOP := -100:
            error "X and Y are not the same size.":
        fi:

        xArr := Array(0..n-1,datatype=integer[8]):
        yArr := Array(0..n-1,datatype=integer[8]):
        for i from 1 to n do
            xArr[i-1] := xVals[i]:
            yArr[i-1] := yVals[i]:
        od:

        nOLEN := N+1:
        dOLEN := DBound+1:
        nOUT := Array(0..nOLEN-1,datatype=integer[8]):
        dOUT := Array(0..dOLEN-1,datatype=integer[8]):

        cppRet := mINTERPDFTRFR(
                          n,xArr,yArr,
                          N,DBound,E,p,
                          nOLEN,nOUT,
                          dOLEN,dOUT
        ):
        lastOP := cppRet:
        if cppRet < 0 then
            error "cppInterpDFTRFR external error code %1", cppRet:
        fi:
        if cppRet = 1 then
            return FAIL:
        fi:

        degN := checkZeroPY(nOUT,nOLEN):
        degD := checkZeroPY(dOUT,dOLEN):
        if degN<0 or degD<0 then return FAIL: fi:

        nn := add(nOUT[i]*var^i, i=0..degN):
        dd := add(dOUT[i]*var^i, i=0..degD):
        return (nn/dd):
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

(* ========================================================================= *)
(* FTRFR bindings : cppDFTRFR (deterministic FT) and cppHFTRFR (gap + bad set) *)
(* REF-FREE: all scalar results come back via the return code and via output  *)
(* arrays (read with checkZeroPY / info[0]), because REF(integer[4]) write-back*)
(* is not reliable here -- the same reason cppRR recomputes degNOUT.           *)
(* Paste ftrfr_cpp_insert.cpp into main.cpp and rebuild cppObj.so first.       *)
(* ========================================================================= *)

mDFTRFR := define_external(
                          'cppDFTRFR',
                          mLen::integer[4],
                          degM::integer[4],
                          M::ARRAY(0..mLen-1, datatype=integer[8]),
                          uLen::integer[4],
                          degU::integer[4],
                          U::ARRAY(0..uLen-1, datatype=integer[8]),
                          N::integer[4],
                          DBound::integer[4],
                          E::integer[4],
                          p::integer[8],
                          nOLEN::integer[4],
                          nOUT::ARRAY(0..nOLEN-1, datatype=integer[8]),
                          dOLEN::integer[4],
                          dOUT::ARRAY(0..dOLEN-1, datatype=integer[8]),
                          RETURN::integer[4],
                          LIB=libObj
                          ):

(* Remove type checking from mDFTRFR. *)
mDFTRFR := subsop(1=(
                    mLen,degM,M,
                    uLen,degU,U,
                    N,DBound,E,p,
                    nOLEN,nOUT,
                    dOLEN,dOUT),
                    op(mDFTRFR)):

mHFTRFR := define_external(
                          'cppHFTRFR',
                          mLen::integer[4],
                          degM::integer[4],
                          M::ARRAY(0..mLen-1, datatype=integer[8]),
                          uLen::integer[4],
                          degU::integer[4],
                          U::ARRAY(0..uLen-1, datatype=integer[8]),
                          p::integer[8],
                          fLEN::integer[4],
                          fOUT::ARRAY(0..fLEN-1, datatype=integer[8]),
                          gLEN::integer[4],
                          gOUT::ARRAY(0..gLEN-1, datatype=integer[8]),
                          lamLEN::integer[4],
                          lamOUT::ARRAY(0..lamLEN-1, datatype=integer[8]),
                          infoLEN::integer[4],
                          infoOUT::ARRAY(0..infoLEN-1, datatype=integer[8]),
                          RETURN::integer[4],
                          LIB=libObj
                          ):

(* Remove type checking from mHFTRFR. *)
mHFTRFR := subsop(1=(
                    mLen,degM,M,
                    uLen,degU,U,
                    p,
                    fLEN,fOUT,
                    gLEN,gOUT,
                    lamLEN,lamOUT,
                    infoLEN,infoOUT),
                    op(mHFTRFR)):


cppDFTRFR := proc(Uin, Min, var, N, DBound, E, p)
    global lastOP:
    local Upoly,Mpoly,degU,degM,uLen,mLen,UArr,MArr,
          nOLEN,dOLEN,nOUT,dOUT,degNOUT,degDOUT,cppRet,nn,dd:

        Upoly := Uin: Mpoly := Min:
        degU := degree(Upoly,var):
        degM := degree(Mpoly,var):
        if degM<0 then
            lastOP := -999:
            error "M is a 0 polynomial.":
        fi:

        mLen := degM+1:
        MArr := convertPY2ARR(Mpoly,var,degM,p):
        if degU<0 then                       (* U == 0 : DFTRFR returns 0/1 *)
            uLen := 1:
            UArr := Array(0..0,datatype=integer[8]):
            degU := -1:
        else
            uLen := degU+1:
            UArr := convertPY2ARR(Upoly,var,degU,p):
        fi:

        nOLEN := N+1:
        dOLEN := DBound+1:
        nOUT := Array(0..nOLEN-1,datatype=integer[8]):
        dOUT := Array(0..dOLEN-1,datatype=integer[8]):

        cppRet := mDFTRFR(
                          mLen,degM,MArr,
                          uLen,degU,UArr,
                          N,DBound,E,p,
                          nOLEN,nOUT,
                          dOLEN,dOUT
        ):
        lastOP := cppRet:

        if cppRet < 0 then
            error "cppDFTRFR external error code %1", cppRet:
        fi:
        if cppRet = 1 then                   (* reconstruction FAILED *)
            return FAIL:
        fi:

        (* degrees recovered from the arrays, not from a REF *)
        degNOUT := checkZeroPY(nOUT,nOLEN):
        degDOUT := checkZeroPY(dOUT,dOLEN):
        if degNOUT<0 or degDOUT<0 then
            return FAIL:
        fi:

        nn := convertARR2PY(nOUT,degNOUT,var):
        dd := convertARR2PY(dOUT,degDOUT,var):
        if dd=0 then
            error "DFTRFR denominator = 0.":
        fi:
        return (nn/dd):
end proc:

cppHFTRFR := proc(Min, Uin, var, p)
    global lastOP:
    local Mpoly,Upoly,degM,degU,mLen,uLen,MArr,UArr,
          fLEN,gLEN,lamLEN,fOUT,gOUT,lamOUT,infoOUT,
          degFOUT,degGOUT,degLamOUT,qmax,cppRet,
          fc,gc,Lambda,rts,badset,i:

        Mpoly := Min: Upoly := Uin:
        degM := degree(Mpoly,var):
        degU := degree(Upoly,var):
        if degM<0 or degU<0 then
            lastOP := -999:
            error "M or U is a 0 polynomial.":
        fi:

        mLen := degM+1: uLen := degU+1:
        MArr := convertPY2ARR(Mpoly,var,degM,p):
        UArr := convertPY2ARR(Upoly,var,degU,p):

        fLEN := degM+1: gLEN := degM+1: lamLEN := degM+1:
        fOUT    := Array(0..fLEN-1,datatype=integer[8]):
        gOUT    := Array(0..gLEN-1,datatype=integer[8]):
        lamOUT  := Array(0..lamLEN-1,datatype=integer[8]):
        infoOUT := Array(0..0,datatype=integer[8]):   (* info[0] = qmax *)

        cppRet := mHFTRFR(
                          mLen,degM,MArr,
                          uLen,degU,UArr,
                          p,
                          fLEN,fOUT,
                          gLEN,gOUT,
                          lamLEN,lamOUT,
                          1,infoOUT
        ):
        lastOP := cppRet:
        if cppRet <> 0 then
            error "cppHFTRFR external error code %1", cppRet:
        fi:

        (* all scalars recovered from arrays -- no REF *)
        degFOUT   := checkZeroPY(fOUT,fLEN):
        degGOUT   := checkZeroPY(gOUT,gLEN):
        degLamOUT := checkZeroPY(lamOUT,lamLEN):
        qmax      := infoOUT[0]:

        fc     := convertARR2PY(fOUT,degFOUT,var):
        gc     := convertARR2PY(gOUT,degGOUT,var):
        Lambda := convertARR2PY(lamOUT,degLamOUT,var):

        if degLamOUT <= 0 then
            badset := []:
        else
            rts := Roots(Lambda) mod p:
            badset := sort([seq(rts[i][1],i=1..nops(rts))]):
        fi:

        return fc,gc,degFOUT,degGOUT,qmax,badset,Lambda:
end proc:

mFTREVAL := define_external(
                          'cppFTREval',
                          nPts::integer[4],
                          alphaLen::integer[4],
                          alpha::ARRAY(0..alphaLen-1, datatype=integer[8]),
                          yLen::integer[4],
                          Yin::ARRAY(0..yLen-1, datatype=integer[8]),
                          sigma::integer[8],
                          N::integer[4],
                          DBound::integer[4],
                          E::integer[4],
                          p::integer[8],
                          outLen::integer[4],
                          out::ARRAY(0..outLen-1, datatype=integer[8]),
                          RETURN::integer[4],
                          LIB=libObj
                          ):

(* Remove type checking from mFTREVAL. *)
mFTREVAL := subsop(1=(nPts,alphaLen,alpha,yLen,Yin,sigma,
                      N,DBound,E,p,outLen,out),
                   op(mFTREVAL)):

cppFTREval := proc(nPts, alphaArr, Yarr, sigma, N, DBound, E, p, outArr)
    local cppRet:
        cppRet := mFTREVAL(nPts,
                           numelems(alphaArr), alphaArr,
                           numelems(Yarr),     Yarr,
                           sigma, N, DBound, E, p,
                           numelems(outArr),   outArr):
        if cppRet < 0 then
            error "cppFTREval external error code %1", cppRet:
        fi:
        if cppRet = 1 then
            return FAIL:
        fi:
        return 0:
end proc:

mAFFINELINE := define_external(
                          'cppAffineLine',
                          T::integer[4],
                          numVar::integer[4],
                          alphaLen::integer[4],
                          alpha::ARRAY(0..alphaLen-1, datatype=integer[8]),
                          betaLen::integer[4],
                          beta::ARRAY(0..betaLen-1, datatype=integer[8]),
                          sigmaLen::integer[4],
                          sigma::ARRAY(0..sigmaLen-1, datatype=integer[8]),
                          p::integer[8],
                          outLen::integer[4],
                          out::ARRAY(0..outLen-1, datatype=integer[8]),
                          RETURN::integer[4],
                          LIB=libObj
                          ):

(* Remove type checking from mAFFINELINE. *)
mAFFINELINE := subsop(1=(T,numVar,alphaLen,alpha,betaLen,beta,
                         sigmaLen,sigma,p,outLen,out),
                      op(mAFFINELINE)):

cppAffineLine := proc(T, numVar, alphaArr, betaArr, sigmaArr, p, outArr)
    global lastOP:
    local cppRet:
        cppRet := mAFFINELINE(T, numVar,
                              numelems(alphaArr), alphaArr,
                              numelems(betaArr),  betaArr,
                              numelems(sigmaArr), sigmaArr,
                              p,
                              numelems(outArr),   outArr):
        lastOP := cppRet:
        if cppRet <> 0 then
            error "cppAffineLine external error code %1", cppRet:
        fi:
        return 0:
end proc:

cppRREFext := define_external("cppRREF",
    n::integer[4],
    m::integer[4],
    B::ARRAY(datatype=integer[8],order=C_order),
    p::integer[8],
    RETURN::integer[4],
    LIB=libObj):

#  Fail here, at read time, rather than 500 black box calls later with
#  "cannot determine if this expression is true or false".
if not type(cppRREFext,procedure) then
    error "define_external did not bind cppRREFext; check that %1 exists and "
          "exports cppRREF (nm -D %1 | grep cppRREF)",POLMATH_LIB:
fi:

#  cppMat64 : build the hardware matrix the external call needs.  Entries are
#  reduced into [0,p) here so that every one of them fits in a machine word.
cppMat64 := proc(A::Matrix,p::prime)
    local n,m,i,j:
    n,m := op(1,A):
    return Matrix(n,m,(i,j) -> modp(A[i,j],p),
                  datatype=integer[8],order=C_order):
end proc:

#  cppRREFip(B,p) puts B in reduced row echelon form IN PLACE and returns the
#  rank.  B must be an integer[8] C_order Matrix.
cppRREFip := proc(B::Matrix(datatype=integer[8]),p::prime)
    local n,m,r:
    n,m := op(1,B):
    r := cppRREFext(n,m,B,p):
    if not type(r,integer) then
        error "cppRREFext returned unevaluated -- the external call never ran; "
              "check LIB=%1 and that pol_math.so exports cppRREF",POLMATH_LIB:
    fi:
    if r < 0 then
        error "cppRREF failed with code %1",r:
    fi:
    return r:
end proc:

#  cppRREF(A,p) returns rref(A) mod p and rank(A).  A itself is left alone.
cppRREF := proc(A::Matrix,p::prime)
    local B,r:
    B := cppMat64(A,p):
    r := cppRREFip(B,p):
    return B,r:
end proc:

#  cppLSip(A,p) solves the n x (n+1) augmented system [A|b] mod p and returns
#  the solution as a list, or FAIL when the system is singular or
#  inconsistent.  A MUST be an integer[8] C_order Matrix and IS DESTROYED --
#  this is the entry point for the black box, where the matrix is built fresh
#  on every call so a defensive copy would be pure overhead.
#
#  After rref the system has a unique solution exactly when the rank is n and
#  the pivots sit in columns 1..n, i.e. A[i,i]=1 for i=1..n.  A pivot in the
#  last column means the system is inconsistent.
cppLSip := proc(A::Matrix(datatype=integer[8]),p::prime)
    local n,m,r,i:
    n,m := op(1,A):
    if m <> n+1 then
        error "expecting an n x (n+1) augmented matrix, got %1 x %2",n,m:
    fi:
    r := cppRREFip(A,p):
    if r <> n then
        return FAIL:
    fi:
    for i from 1 to n do
        if A[i,i] <> 1 then
            return FAIL:
        fi:
    od:
    return [seq(A[i,n+1],i=1..n)]:
end proc:

#  cppLS(A,p) is the same thing for a general Matrix.  It copies first, so the
#  caller's matrix survives the call.
cppLS := proc(A::Matrix,p::prime)
    return cppLSip(cppMat64(A,p),p):
end proc: