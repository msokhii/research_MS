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

cppEvalSolveext := define_external("cppEvalSolve",
    nr::integer[4],
    nc::integer[4],
    nv::integer[4],
    entStart::ARRAY(datatype=integer[4]),
    expo::ARRAY(datatype=integer[4]),
    coef::ARRAY(datatype=integer[8]),
    pnt::ARRAY(datatype=integer[8]),
    p::integer[8],
    dmax::integer[4],
    outLen::integer[4],
    xOut::ARRAY(datatype=integer[8]),
    infoLen::integer[4],
    info::ARRAY(datatype=integer[8]),
    RETURN::integer[4],
    LIB=libObj):

cppEvalSolveBlockext := define_external("cppEvalSolveBlock",
    nr::integer[4],
    nc::integer[4],
    nv::integer[4],
    entStart::ARRAY(datatype=integer[4]),
    expo::ARRAY(datatype=integer[4]),
    coef::ARRAY(datatype=integer[8]),
    pts::ARRAY(datatype=integer[8]),
    npts::integer[4],
    p::integer[8],
    dmax::integer[4],
    ldx::integer[4],
    xOut::ARRAY(datatype=integer[8]),
    infoLen::integer[4],
    info::ARRAY(datatype=integer[8]),
    RETURN::integer[4],
    LIB=libObj):

#cppRootsext := define_external("cppRoots",
#    degF::integer[4],
#    f::ARRAY(datatype=integer[8]),
#    p::integer[8],
#    outLen::integer[4],
#    rootsOut::ARRAY(datatype=integer[8]),
#    multLen::integer[4],
#    multOut::ARRAY(datatype=integer[8]),
#    infoLen::integer[4],
#    info::ARRAY(datatype=integer[8]),
#    RETURN::integer[4],
#    LIB=POLMATH_LIB):

#  Fail here, at read time, rather than 500 black box calls later with
#  "cannot determine if this expression is true or false".
if not type(cppRREFext,procedure) then
    error "define_external did not bind cppRREFext; check that %1 exists and "
          "exports cppRREF (nm -D %1 | grep cppRREF)",POLMATH_LIB:
fi:
if not type(cppEvalSolveext,procedure) then
    error "define_external did not bind cppEvalSolveext; check that %1 exists "
          "and exports cppEvalSolve",POLMATH_LIB:
fi:
if not type(cppEvalSolveBlockext,procedure) then
    error "define_external did not bind cppEvalSolveBlockext; check that %1 "
          "exists and exports cppEvalSolveBlock",POLMATH_LIB:
fi:
#if not type(cppRootsext,procedure) then
#    error "define_external did not bind cppRootsext; check that %1 exists and "
#          "exports cppRoots",POLMATH_LIB:
#fi:

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
    #  If B were Fortran order, define_external would silently hand the C code
    #  a transposed COPY, the rref would land in that copy, and B would come
    #  back untouched.  Refuse instead of returning nonsense.
    if rtable_options(B,order) <> C_order then
        error "expecting a C_order Matrix":
    fi:
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

#  ---------------------------------------------------------------------------
#  EVALUATING THE PARAMETRIC MATRIX
#
#  modp(eval(L[i,j],sv),p) evaluates over Z FIRST and only then reduces.  With
#  64 bit points and entries of degree d in the parameters the intermediates
#  are ~62*d bit integers, so the cost is bignum arithmetic, not the solve.
#  Both routes below reduce as they go.
#
#  Only Modular:-LinearSolve is capped at 2^32; Modular:-Mod itself may well
#  accept a 64 bit modulus with datatype=integer[8], in which case the
#  evaluation runs in exactly the compiled kernel Maple was using before.
#  The first call finds out and the answer is cached in CPPEVALMODE.
#     CPPEVALMODE = 0   not yet decided
#                 = 1   LinearAlgebra:-Modular:-Mod
#                 = 2   Eval ... mod p
#  ---------------------------------------------------------------------------

#  Maple's modp1 / Modular kernels only cover primes below about 2^31.5, and
#  below that bound they beat the external route.  Above it they fall back to
#  generic Zp arithmetic and the C++ path wins.  Every dispatch in this file
#  uses this one constant.
MAPLE_FAST_BOUND := 2^31:

CPPEVALMODE := 0:

cppEvalMatrix64 := proc(L::Matrix,sv::list,p::prime)
    local A,nr,nc,svs,i,j:
    global CPPEVALMODE:
    nr,nc := op(1,L):
    if CPPEVALMODE = 0 then
        A := traperror(LinearAlgebra:-Modular:-Mod(p,L,sv,integer[8])):
        if type(A,Matrix) and rtable_options(A,datatype) = integer[8]
                          and rtable_options(A,order) = C_order then
            CPPEVALMODE := 1:
            return A:
        fi:
        CPPEVALMODE := 2:
        WARNING("Modular:-Mod would not build an integer[8] C_order matrix "
                "for this modulus; falling back to Eval ... mod p"):
    fi:
    if CPPEVALMODE = 1 then
        return LinearAlgebra:-Modular:-Mod(p,L,sv,integer[8]):
    fi:
    svs := {op(sv)}:
    return Matrix(nr,nc,(i,j) -> Eval(L[i,j],svs) mod p,
                  datatype=integer[8],order=C_order):
end proc:

#  ---------------------------------------------------------------------------
#  PARAMETRIC MATRIX EVALUATION + SOLVE
#
#  The black box evaluates ONE fixed augmented matrix L(y1,..,ynv) at many
#  points.  Doing that in Maple costs an interpreted eval plus a modp per
#  entry, which measured at about 1.5 microseconds per entry -- more than the
#  linear solve itself.  cppEncodeMatrix turns L into a flat sparse monomial
#  encoding ONCE, and cppEvalSolve then does evaluation and solve in a single
#  external call.  Per call Maple only fills the nv parameter values.
#
#  The encoding is CSR over the nr*nc entries in row major order:
#    ent[e] .. ent[e+1]-1   are the term indices of entry e
#    coefRaw[t+1]           exact Maple coefficient of term t
#    coef[t]                the same reduced mod p (rebuilt when p changes)
#    expo[t*nv+k]           exponent of params[k+1] in term t
#  ---------------------------------------------------------------------------

cppEncodeMatrix := proc(L::Matrix,params::list)
    local nr,nc,nv,i,j,k,t,e,cf,mm,T,nT,dmax,d,degs,idx,ent,expoA,coefL,E:
    nr,nc := op(1,L):
    nv := nops(params):
    if nv < 1 then
        error "expecting at least one parameter":
    fi:
    T    := table():
    nT   := 0:
    dmax := 0:
    for i from 1 to nr do
        for j from 1 to nc do
            e := expand(L[i,j]):
            if not type(e,polynom(rational,params)) then
                error "entry (%1,%2) must be a polynomial in the parameters "
                      "with rational coefficients",i,j:
            fi:
            if e = 0 then
                T[i,j] := [[],[]]:
                next:
            fi:
            cf := [coeffs(e,params,'mm')]:
            mm := [mm]:
            degs := []:
            for k from 1 to nops(cf) do
                for t from 1 to nv do
                    d := degree(mm[k],params[t]):
                    if d > dmax then dmax := d: fi:
                    degs := [op(degs),d]:
                od:
            od:
            T[i,j] := [cf,degs]:
            nT := nT+nops(cf):
        od:
    od:
    ent   := Array(0..nr*nc,datatype=integer[4]):
    expoA := Array(0..max(nv*nT,1)-1,datatype=integer[4]):
    coefL := Array(1..max(nT,1)):                 # exact, may be big or rational
    idx := 0:
    for i from 1 to nr do
        for j from 1 to nc do
            ent[(i-1)*nc+j-1] := idx:
            cf   := T[i,j][1]:
            degs := T[i,j][2]:
            for k from 1 to nops(cf) do
                coefL[idx+1] := cf[k]:
                for t from 1 to nv do
                    expoA[idx*nv+t-1] := degs[(k-1)*nv+t]:
                od:
                idx := idx+1:
            od:
        od:
    od:
    ent[nr*nc] := nT:

    #  Everything the external call needs, allocated once.  The work arrays are
    #  reused on every black box call so the hot path does not allocate.
    E := table():
    E["L"]       := L:
    E["params"]  := params:
    E["nr"]      := nr:
    E["nc"]      := nc:
    E["nv"]      := nv:
    E["nT"]      := nT:
    E["dmax"]    := dmax:
    E["ent"]     := ent:
    E["expo"]    := expoA:
    E["coefRaw"] := coefL:
    E["p"]       := 0:
    E["coef"]    := Array(0..max(nT,1)-1,datatype=integer[8]):
    E["pnt"]     := Array(0..nv-1,datatype=integer[8]):
    E["x"]       := Array(0..nr-1,datatype=integer[8]):
    E["info"]    := Array(0..1,datatype=integer[8]):
    return E:
end proc:

#  cppEvalSolve(E,point_,p) evaluates the encoded matrix at point_ and solves
#  the system mod p.  Returns the solution as a list, or FAIL when the
#  evaluated matrix is singular or the system is inconsistent.
#  The coefficients are reduced mod p only when p changes.
cppEvalSolve := proc(E,point_::list,p::prime)
    local k,rc,nr,nv,A,T:

    #  Below the bound Maple's own compiled path wins: measured at n=12 it is
    #  2.2x faster per black box call than going out to C++.  Only use the
    #  external route where Maple cannot follow.
    if p < MAPLE_FAST_BOUND then
        A := LinearAlgebra:-Modular:-Mod(p,E["L"],
                 zip((par,pnt) -> par=pnt,E["params"],point_),integer):
        T := traperror(LinearAlgebra:-Modular:-LinearSolve(p,A,1)):
        if T = "Matrix is singular." then
            return FAIL:
        fi:
        return convert(A[1..E["nr"],E["nc"]],list):
    fi:

    nr := E["nr"]:
    nv := E["nv"]:
    if nops(point_) < nv then
        error "expecting %1 parameter values, got %2",nv,nops(point_):
    fi:
    if E["p"] <> p then
        for k from 1 to E["nT"] do
            E["coef"][k-1] := modp(E["coefRaw"][k],p):
        od:
        E["p"] := p:
    fi:
    #  One kernel conversion instead of nv interpreted modp + element stores.
    #  The C side reduces anything outside [0,p) itself.
    E["pnt"] := Array(point_,datatype=integer[8]):
    rc := cppEvalSolveext(nr,E["nc"],nv,E["ent"],E["expo"],E["coef"],
                          E["pnt"],p,E["dmax"],nr,E["x"],2,E["info"]):
    if not type(rc,integer) then
        error "cppEvalSolveext returned unevaluated -- the external call never "
              "ran; check LIB=%1",POLMATH_LIB:
    fi:
    if rc < 0 then
        error "cppEvalSolve failed with code %1",rc:
    fi:
    if rc = 1 then
        return FAIL:
    fi:
    return convert(E["x"],list):
end proc:

#  ---------------------------------------------------------------------------
#  ROOTS OVER GF(p)   -- COMMENTED OUT, Roots(F) mod p is used instead.
#
#  At n=12 the inputs are degree 2563 and 5579 and the classical O(d^2 log p)
#  Cantor-Zassenhaus below loses to Maple.  cppRoots stays in pol_math.cpp,
#  unreferenced; uncomment this section and the three call sites in
#  solver_pol_sys.mpl to put it back.
#
#  Drop in replacement for  Roots(F) mod p .  Same shape of answer: a list of
#  [root,multiplicity] pairs, sorted by root, so R[1][1] = 0 still detects the
#  zero root.  Maple's Roots falls off its fast univariate representation above
#  p = 2^31.5, which is why it costs 14x more at a 64 bit prime than a 32 bit
#  one; the C++ side is built on polmul64s and polDIVIP64, which measured the
#  same at both sizes.
#  ---------------------------------------------------------------------------

#cppRootsOf := proc(F::polynom,x::name,p::prime)
#    local d,fA,rts,mlt,info,rc,n,k:
#    d := degree(F,x):
#    if d < 1 then
#        return []:
#    fi:
#    #  Same dispatch as cppEvalSolve.  Measured at n=12 with degree 2563 and
#    #  5579 inputs, Maple's Roots is 2.4x faster than the C++ Cantor-Zassenhaus
#    #  below the bound and 13.7x slower above it.
#    if p < MAPLE_FAST_BOUND then
#        return Roots(F) mod p:
#    fi:
#    fA  := Array(0..d,[seq(modp(coeff(F,x,k),p),k=0..d)],datatype=integer[8]):
#    rts := Array(0..d-1,datatype=integer[8]):
#    mlt := Array(0..d-1,datatype=integer[8]):
#    info := Array(0..0,datatype=integer[8]):
#    rc := cppRootsext(d,fA,p,d,rts,d,mlt,1,info):
#    if not type(rc,integer) then
#        error "cppRootsext returned unevaluated -- the external call never "
#              "ran; check LIB=%1",POLMATH_LIB:
#    fi:
#    if rc < 0 then
#        error "cppRoots failed with code %1",rc:
#    fi:
#    n := info[0]:
#    return [seq([rts[k],mlt[k]],k=0..n-1)]:
#end proc:

#  ---------------------------------------------------------------------------
#  BATCHED, TRANSPOSED BLACK BOX
#
#  The MRFI loop wants, for each equation i, the sequence of x_i over the whole
#  block of points on the affine line.  Point at a time that costs npts
#  external calls, npts solution lists, and an interpreted transpose:
#
#      BBvals := [seq(B(Psi_alpha[s],p),s=1..mMax)]:
#      for i ... for s ... Yarr[s-1] := BBvals[s][i]
#
#  cppEvalSolveBlock does the whole block in ONE call and writes the answers
#  already transposed, so row i is contiguous and ArrayTools:-Alias hands it to
#  cppFTREval with no copy.  Nothing is built and then rearranged.
#
#  E   from cppEncodeMatrix
#  pts flat Array of npts*nv values, point major, exactly as cppAffineLine
#      writes them (get_point_block_on_affine_line returns this)
#
#  Returns an nr x npts row major Array X with X[(i-1)*npts+(s-1)] = x_i at
#  point s, or FAIL when some point gave a singular system.
#  ---------------------------------------------------------------------------

cppEvalSolveBlock := proc(E,pts::Array,npts::posint,p::prime)
    local k,s,rc,nr,nv,X,soln:
    nr := E["nr"]:
    nv := E["nv"]:
    if E["p"] <> p then
        for k from 1 to E["nT"] do
            E["coef"][k-1] := modp(E["coefRaw"][k],p):
        od:
        E["p"] := p:
    fi:
    #  One buffer per (nr,npts) shape, reused across calls.
    if not assigned(E["Xn"]) or E["Xn"] < nr*npts then
        E["X"]  := Array(0..nr*npts-1,datatype=integer[8]):
        E["Xn"] := nr*npts:
    fi:
    X := E["X"]:

    #  Below the bound cppEvalSolve takes Maple's compiled Modular path, which
    #  is faster per solve than going out to C++.  Loop it, but still write the
    #  answers transposed so the caller gets the same contiguous rows.
    if p < MAPLE_FAST_BOUND then
        for s from 1 to npts do
            soln := cppEvalSolve(E,[seq(pts[(s-1)*nv+k],k=0..nv-1)],p):
            if soln = FAIL then
                E["info"][1] := s:
                return FAIL:
            fi:
            for k from 1 to nr do
                X[(k-1)*npts+s-1] := soln[k]:
            od:
        od:
        return X:
    fi:

    rc := cppEvalSolveBlockext(nr,E["nc"],nv,E["ent"],E["expo"],E["coef"],
                               pts,npts,p,E["dmax"],npts,X,2,E["info"]):
    if not type(rc,integer) then
        error "cppEvalSolveBlockext returned unevaluated -- the external call "
              "never ran; check LIB=%1",POLMATH_LIB:
    fi:
    if rc < 0 then
        error "cppEvalSolveBlock failed with code %1",rc:
    fi:
    if rc = 1 then
        return FAIL:                      # E["info"][1] is the first bad point
    fi:
    return X:
end proc:

#  Row i of the block as a zero copy alias.  m entries starting at x_i(point 1).
#  The result shares storage with X, so writing into it (the fault injection
#  path does) writes into X -- that is harmless here because row i is finished
#  with before the next equation is touched, but it is worth knowing.
cppBlockRow := proc(X::Array,i::posint,npts::posint,m::posint)
    #  0..m-1 on purpose: cppFTREval and the fault path index Yarr[s-1].
    return ArrayTools:-Alias(X,(i-1)*npts,[0..m-1]):
end proc:

#  ---------------------------------------------------------------------------
#  MONAGAN'S polroots64s THROUGH THE FFI
#
#  cppPolRoots is now part of pol_math.cpp, so this is on.  Set it back to false
#  to unbind everything here and fall back to Roots(F) mod p without touching
#  the call sites.
#
#  ROOTS_MODE picks the dispatch:
#     0  by prime size: Maple below MAPLE_FAST_BOUND, polroots64s above
#     1  always Roots(F) mod p
#     2  always polroots64s
#  ---------------------------------------------------------------------------

USE_POLROOTS64S := true:
ROOTS_MODE      := 0:

if USE_POLROOTS64S then
    cppPolRootsext := define_external("cppPolRoots",
        d::integer[4],
        f::ARRAY(datatype=integer[8]),
        p::integer[8],
        wsize::integer[4],
        outLen::integer[4],
        rootsOut::ARRAY(datatype=integer[8]),
        infoLen::integer[4],
        info::ARRAY(datatype=integer[8]),
        RETURN::integer[4],
        LIB=libObj):
    if not type(cppPolRootsext,procedure) then
        error "define_external did not bind cppPolRootsext; check that %1 "
              "exports cppPolRoots",POLMATH_LIB:
    fi:
fi:

#  polroots64s returns DISTINCT roots, so every multiplicity here is 1.  The
#  shape matches Roots(F) mod p: a list of [root,multiplicity] pairs sorted by
#  root, so R[1][1] = 0 still detects the zero root.
cppPolRootsOf := proc(F::polynom,x::name,p::prime,wsize::integer := 0)
    local d,fA,rts,info,rc,n,k:
    d := degree(F,x):
    if d < 1 then
        return []:
    fi:
    fA   := Array(0..d,[seq(modp(coeff(F,x,k),p),k=0..d)],datatype=integer[8]):
    rts  := Array(0..d-1,datatype=integer[8]):
    info := Array(0..0,datatype=integer[8]):
    rc := cppPolRootsext(d,fA,p,wsize,d,rts,1,info):
    if not type(rc,integer) then
        error "cppPolRootsext returned unevaluated -- the external call never "
              "ran; check LIB=%1",POLMATH_LIB:
    fi:
    if rc = 1 then
        error "cppPolRoots found %1 roots, more than the %2 slots given",
              info[0],d:
    fi:
    if rc < 0 then
        error "cppPolRoots failed with code %1",rc:
    fi:
    n := info[0]:
    return [seq([rts[k],1],k=0..n-1)]:
end proc:

#  Single dispatch point for the three MRFI call sites.
rootsMODp := proc(F::polynom,x::name,p::prime)
    if ROOTS_MODE = 1 or not USE_POLROOTS64S then
        return Roots(F) mod p:
    fi:
    if ROOTS_MODE = 2 or p >= MAPLE_FAST_BOUND then
        return cppPolRootsOf(F,x,p):
    fi:
    return Roots(F) mod p:
end proc:
