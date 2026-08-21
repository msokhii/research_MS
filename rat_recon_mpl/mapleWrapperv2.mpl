libObj := "/cecm/home/mss59/Desktop/resMaple_MS/inuse_cpp_routines/cpp_routines/cppObj_mserver.so":
#libObj := "/home/msokhi/Desktop/res_MS/inuse_cpp_routines/cpp_routines/cppObj.so":

# =============================================================================
# External procedure: mRATRECON
# What it does:
#   Binds Maple directly to the C++ ratRECON_C rational-reconstruction entry point.
# Inputs:
#   - mLen,degM,M: modulus-polynomial length, degree, and integer[8] coefficients.
#   - uLen,degU,U: interpolant length, degree, and coefficients.
#   - N,DBound,p: degree bounds and modulus.
#   - nOUT,dOUT plus lengths and degree references: caller-provided outputs.
# Outputs:
#   - Returns the C++ integer status code and writes numerator/denominator coefficients and degrees into the supplied outputs.
# Example:
#   rc := mRATRECON(mLen,degM,M,uLen,degU,U,N,D,p,nLen,nOut,dn,dLen,dOut,dd):
# =============================================================================
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

# =============================================================================
# External procedure: mNEWTONINTERP
# What it does:
#   Binds Maple directly to the C++ cppInterp Newton interpolation wrapper.
# Inputs:
#   - xLen,xIn and yLen,yIn: integer[8] sample arrays.
#   - p: modulus.
#   - outLen,yOut: coefficient output buffer.
#   - degOut: output degree reference.
# Outputs:
#   - Returns the C++ status code and fills yOut/degOut on success.
# Example:
#   rc := mNEWTONINTERP(n,xA,n,yA,p,n,cA,d):
# =============================================================================
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

# =============================================================================
# External procedure: mVSOLVE
# What it does:
#   Binds Maple directly to the C++ transposed Vandermonde solver.
# Inputs:
#   - mLen,mIn: node array.
#   - yLen,yIn: right-hand side.
#   - shiftInt: shift correction.
#   - pp: modulus.
#   - outLen,aOut: solution buffer.
# Outputs:
#   - Returns the C++ status code and writes the coefficient solution into aOut.
# Example:
#   rc := mVSOLVE(n,mA,n,yA,0,p,n,aA):
# =============================================================================
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

# =============================================================================
# External procedure: mBM
# What it does:
#   Binds Maple directly to the C++ Berlekamp-Massey wrapper.
# Inputs:
#   - aLen,aIn: input sequence.
#   - p: modulus.
#   - outLen,lOut: connection-polynomial buffer.
#   - degOut: output degree reference.
# Outputs:
#   - Returns the C++ status code and fills lOut/degOut.
# Example:
#   rc := mBM(N,aA,p,n+1,L,d):
# =============================================================================
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

# =============================================================================
# Procedure: cppBMM
# What it does:
#   Wraps the C++ Berlekamp-Massey routine and converts its output coefficient array into a Maple list.
# Inputs:
#   - a: nonempty sequence supplied as a Vector or list.
#   - p: prime modulus.
# Outputs:
#   - Returns the connection-polynomial coefficients [L[0],...,L[d]] in ascending degree order; returns [] when the C++ routine produces no polynomial.
# Example:
#   LambdaCoeffs := cppBMM([1,2,4,8],101):
# =============================================================================
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

# =============================================================================
# Procedure: cppVS
# What it does:
#   Wraps the C++ transposed Vandermonde solver used by sparse interpolation.
# Inputs:
#   - v: right-hand-side/sample values.
#   - m: Vandermonde nodes; must have the same length as v.
#   - p: prime modulus.
#   - shift: optional integer shift correction, default 0.
# Outputs:
#   - Returns the solved coefficient list of the same length as v.
# Example:
#   a := cppVS([3,5],[2,7],101):
# =============================================================================
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

# =============================================================================
# External procedure: mINTERPDFTRFR
# What it does:
#   Binds Maple to the combined C++ interpolation plus deterministic fault-tolerant reconstruction routine.
# Inputs:
#   - nPts,alpha,Yin: sample count, nodes, and values.
#   - N,DBound,E,p: degree/error bounds and modulus.
#   - nOUT,dOUT with lengths: output coefficient arrays.
# Outputs:
#   - Returns 0 on successful reconstruction, 1 on reconstruction failure, or a negative C++ error code.
# Example:
#   rc := mINTERPDFTRFR(n,alphaA,yA,N,D,E,p,N+1,nA,D+1,dA):
# =============================================================================
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


# =============================================================================
# Procedure: cppInterpDFTRFR
# What it does:
#   Performs univariate Newton interpolation followed immediately by deterministic fault-tolerant rational reconstruction through the C++ combined wrapper.
# Inputs:
#   - xVals,yVals: equally sized sample-node and sample-value lists.
#   - var: Maple name used to build the returned rational function.
#   - N,DBound: numerator and denominator degree bounds.
#   - E: error budget.
#   - p: modulus.
# Outputs:
#   - Returns the reconstructed rational expression in var, or FAIL when reconstruction does not validate.
# Example:
#   f := cppInterpDFTRFR([1,2,3],[4,7,12],x,1,1,0,101):
# =============================================================================
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
# =============================================================================
# Procedure: checkZeroPY
# What it does:
#   Finds the highest nonzero coefficient index in a raw coefficient Array.
# Inputs:
#   - A: coefficient Array indexed from 0.
#   - len: number of slots to inspect.
# Outputs:
#   - Returns the largest i in 0..len-1 with A[i]<>0, or -1 for the zero polynomial.
# Example:
#   d := checkZeroPY(A,numelems(A)):
# =============================================================================
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
# =============================================================================
# Procedure: convertPY2ARR
# What it does:
#   Copies the coefficients of a Maple polynomial into an integer[8] Array ordered from constant term to leading term.
# Inputs:
#   - poly: polynomial expression.
#   - var: polynomial variable.
#   - deg: degree to copy.
#   - p: modulus argument kept for interface consistency; this procedure does not reduce coefficients itself.
# Outputs:
#   - Returns Array(0..deg) with A[i]=coeff(poly,var,i).
# Example:
#   A := convertPY2ARR(3+2*x+x^2,x,2,101):
# =============================================================================
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
# =============================================================================
# Procedure: convertARR2PY
# What it does:
#   Reconstructs a Maple polynomial from an ascending coefficient Array.
# Inputs:
#   - A: coefficient Array indexed from 0.
#   - deg: degree to include.
#   - var: polynomial variable.
# Outputs:
#   - Returns add(A[i]*var^i,i=0..deg).
# Example:
#   f := convertARR2PY(Array(0..2,[3,2,1]),2,x):
# =============================================================================
convertARR2PY := proc(A,deg,var) option inline:
    local i:
    
    add(A[i]*var^i,i=0..deg):
end proc:

(* Check if deg(A[i])=0 *)
# =============================================================================
# Procedure: checkZeroPY
# What it does:
#   Finds the highest nonzero coefficient index in a raw coefficient Array; this later definition also prints a diagnostic for the zero polynomial.
# Inputs:
#   - A: coefficient Array indexed from 0.
#   - len: number of slots to inspect.
# Outputs:
#   - Returns the largest nonzero index, or -1 after printing "ZERO POLYNOMIAL." when all entries vanish.
# Example:
#   d := checkZeroPY(A,numelems(A)):
# =============================================================================
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
# =============================================================================
# Procedure: cppRR
# What it does:
#   Maple-facing rational-reconstruction wrapper around mRATRECON. It converts U and M to coefficient arrays, calls C++, and converts the normalized result back to Maple form.
# Inputs:
#   - Uin,Min: univariate polynomials satisfying the reconstruction congruence.
#   - var: polynomial variable.
#   - N,DBound: numerator and denominator degree bounds.
#   - p: modulus.
# Outputs:
#   - Returns the reconstructed rational function, or FAIL when the C++ reconstruction returns a nonzero status.
# Example:
#   r := cppRR(U,M,x,degN,degD,101):
# =============================================================================
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
# =============================================================================
# Procedure: cppNewtonInterp
# What it does:
#   Interpolates a univariate polynomial from sample points using the C++ Newton interpolation kernel.
# Inputs:
#   - xVals,yVals: equally sized nonempty lists of nodes and values.
#   - var: variable for the returned polynomial.
#   - p: modulus.
# Outputs:
#   - Returns the interpolating polynomial modulo p, or FAIL if the C++ call fails.
# Example:
#   f := cppNewtonInterp([0,1,2],[1,2,5],x,101):  # 1+x^2
# =============================================================================
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

# =============================================================================
# External procedure: mDFTRFR
# What it does:
#   Binds Maple directly to the C++ deterministic fault-tolerant rational reconstruction wrapper.
# Inputs:
#   - M/U lengths, degrees, and coefficient Arrays.
#   - N,DBound,E,p: reconstruction bounds and modulus.
#   - nOUT,dOUT: output buffers with lengths.
# Outputs:
#   - Returns the C++ status code and writes numerator/denominator coefficients into the output Arrays.
# Example:
#   rc := mDFTRFR(mLen,dM,M,uLen,dU,U,N,D,E,p,N+1,nA,D+1,dA):
# =============================================================================
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

# =============================================================================
# External procedure: mHFTRFR
# What it does:
#   Binds Maple directly to the C++ heuristic/gap fault-tolerant reconstruction wrapper.
# Inputs:
#   - M/U arrays with lengths/degrees and p.
#   - fOUT,gOUT,lamOUT: candidate and bad-factor output Arrays.
#   - infoOUT: metadata output, with infoOUT[0]=qmax.
# Outputs:
#   - Returns the C++ status code and fills all supplied output Arrays.
# Example:
#   rc := mHFTRFR(mLen,dM,M,uLen,dU,U,p,fLen,fA,gLen,gA,lLen,lA,1,info):
# =============================================================================
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


# =============================================================================
# Procedure: cppDFTRFR
# What it does:
#   Maple wrapper for deterministic fault-tolerant rational reconstruction from U modulo M.
# Inputs:
#   - Uin,Min: univariate polynomials.
#   - var: polynomial variable.
#   - N,DBound: numerator/denominator degree bounds.
#   - E: error budget.
#   - p: modulus.
# Outputs:
#   - Returns the normalized reconstructed rational function; returns FAIL when the deterministic degree/gcd checks fail.
# Example:
#   r := cppDFTRFR(U,M,x,N,D,E,101):
# =============================================================================
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

# =============================================================================
# Procedure: cppHFTRFR
# What it does:
#   Maple wrapper for the heuristic/gap phase of fault-tolerant rational reconstruction and bad-point detection.
# Inputs:
#   - Min,Uin: nonzero univariate polynomials.
#   - var: polynomial variable.
#   - p: modulus.
# Outputs:
#   - Returns fc,gc,degF,degG,qmax,badset,Lambda, where Lambda is the common bad factor and badset contains its roots modulo p.
# Example:
#   fc,gc,dF,dG,qmax,bad,L := cppHFTRFR(M,U,x,101):
# =============================================================================
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

# =============================================================================
# External procedure: mFTREVAL
# What it does:
#   Binds Maple to the C++ routine that interpolates, reconstructs, and evaluates numerator/denominator values at sigma.
# Inputs:
#   - nPts,alpha,Yin: active sample data.
#   - sigma,N,DBound,E,p: evaluation/reconstruction parameters.
#   - outLen,out: output Array.
# Outputs:
#   - Returns the C++ status code; on success out[0] and out[1] hold numerator(sigma) and denominator(sigma).
# Example:
#   rc := mFTREVAL(n,n,alphaA,n,yA,sigma,N,D,E,p,2,outA):
# =============================================================================
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

# =============================================================================
# Procedure: cppFTREval
# What it does:
#   Calls the combined C++ interpolation/reconstruction/evaluation kernel for one sample sequence.
# Inputs:
#   - nPts: number of active points.
#   - alphaArr,Yarr: node and value Arrays.
#   - sigma: point where the reconstructed numerator/denominator are evaluated.
#   - N,DBound,E,p: reconstruction bounds and modulus.
#   - outArr: output Array supplied by the caller.
# Outputs:
#   - Returns 0 on success or FAIL on reconstruction failure; the C++ call writes numerator(sigma) and denominator(sigma) into outArr.
# Example:
#   status := cppFTREval(n,alphaArr,Yarr,sigma,N,D,E,p,outArr):
# =============================================================================
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

# =============================================================================
# External procedure: mAFFINELINE
# What it does:
#   Binds Maple to the C++ affine-line point generator.
# Inputs:
#   - T,numVar: point and parameter counts.
#   - alpha,beta,sigma Arrays with lengths.
#   - p: modulus.
#   - outLen,out: flat point-major output buffer.
# Outputs:
#   - Returns the C++ status code and fills out with T*numVar residues.
# Example:
#   rc := mAFFINELINE(T,nv,T,aA,nv-1,bA,nv,sA,p,T*nv,outA):
# =============================================================================
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

# =============================================================================
# Procedure: cppAffineLine
# What it does:
#   Calls the C++ affine-line generator and fills a preallocated flat point-major Array.
# Inputs:
#   - T: number of points.
#   - numVar: number of parameters.
#   - alphaArr,betaArr,sigmaArr: line-definition Arrays.
#   - p: modulus.
#   - outArr: output Array of length at least T*numVar.
# Outputs:
#   - Returns 0 on success; outArr receives all generated parameter points.
# Example:
#   cppAffineLine(T,nv,alphaA,betaA,sigmaA,p,ptsA):
# =============================================================================
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

# =============================================================================
# External procedure: cppRREFext
# What it does:
#   Raw external binding to the C++ in-place RREF routine.
# Inputs:
#   - n,m: matrix dimensions.
#   - B: C_order integer[8] Matrix.
#   - p: modulus.
# Outputs:
#   - Returns rank(B) or a negative error code; B is overwritten with its RREF.
# Example:
#   r := cppRREFext(n,m,B,p):
# =============================================================================
cppRREFext := define_external("cppRREF",
    n::integer[4],
    m::integer[4],
    B::ARRAY(datatype=integer[8],order=C_order),
    p::integer[8],
    RETURN::integer[4],
    LIB=libObj):

cppRREFext := subsop(1=(
                       n,
                       m,
                       B,
                       p),
                       op(cppRREFext)):

# =============================================================================
# External procedure: cppEvalSolveext
# What it does:
#   Raw external binding that evaluates one encoded parametric augmented matrix and solves it.
# Inputs:
#   - Dimensions and sparse encoding Arrays entStart/expo/coef.
#   - pnt: parameter point.
#   - p,dmax: modulus and maximum exponent.
#   - xOut and info: solution and rank outputs.
# Outputs:
#   - Returns 0 for a unique solution, 1 for singular/inconsistent input, or a negative error code.
# Example:
#   rc := cppEvalSolveext(nr,nc,nv,ent,expo,coef,pnt,p,dmax,nr,x,2,info):
# =============================================================================
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

cppEvalSolveext := subsop(1=(
                       nr,
                       nc,
                       nv,
                       entStart,
                       expo,
                       coef,
                       pnt,
                       p,
                       dmax,
                       outLen,
                       xOut,
                       infoLen,
                       info),
                       op(cppEvalSolveext)):

# =============================================================================
# External procedure: cppEvalSolveBlockext
# What it does:
#   Raw external binding for batched encoded-matrix evaluation and solving.
# Inputs:
#   - Dimensions and sparse encoding Arrays.
#   - pts,npts: flat block of parameter points.
#   - p,dmax,ldx: arithmetic and output-stride settings.
#   - xOut,info: transposed solutions and status metadata.
# Outputs:
#   - Returns 0 if all points solve, 1 if any point fails, or a negative error code.
# Example:
#   rc := cppEvalSolveBlockext(nr,nc,nv,ent,expo,coef,pts,T,p,dmax,T,X,2,info):
# =============================================================================
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

cppEvalSolveBlockext := subsop(1=(
                       nr,
                       nc,
                       nv,
                       entStart,
                       expo,
                       coef,
                       pts,
                       npts,
                       p,
                       dmax,
                       ldx,
                       xOut,
                       infoLen,
                       info),
                       op(cppEvalSolveBlockext)):

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
#  cppMat64 : build the hardware matrix the external call needs.  Entries are
#  reduced into [0,p) here so that every one of them fits in a machine word.
# =============================================================================
# Procedure: cppMat64
# What it does:
#   Converts a general Maple Matrix into a C-order integer[8] Matrix with every entry reduced modulo p.
# Inputs:
#   - A: source Matrix.
#   - p: prime modulus.
# Outputs:
#   - Returns a new Matrix suitable for the C++ RREF/linear-solve wrappers; A is unchanged.
# Example:
#   B := cppMat64(A,101):
# =============================================================================
cppMat64 := proc(A::Matrix,p::prime)
    local n,m,i,j:
    n,m := op(1,A):
    return Matrix(n,m,(i,j) -> modp(A[i,j],p),
                  datatype=integer[8],order=C_order):
end proc:

#  cppRREFip(B,p) puts B in reduced row echelon form IN PLACE and returns the
#  rank.  B must be an integer[8] C_order Matrix.
# =============================================================================
# Procedure: cppRREFip
# What it does:
#   Runs the C++ reduced-row-echelon-form routine directly on a C-order integer[8] Matrix.
# Inputs:
#   - B: Matrix(datatype=integer[8]) in C_order; modified in place.
#   - p: prime modulus.
# Outputs:
#   - Returns rank(B); B is replaced by its RREF modulo p.
# Example:
#   r := cppRREFip(B,101):
# =============================================================================
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
# =============================================================================
# Procedure: cppRREF
# What it does:
#   Convenience wrapper that copies/converts a general Matrix and computes its modular RREF without changing the caller's matrix.
# Inputs:
#   - A: source Matrix.
#   - p: prime modulus.
# Outputs:
#   - Returns B,r where B is rref(A) mod p and r is its rank.
# Example:
#   B,r := cppRREF(A,101):
# =============================================================================
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
# =============================================================================
# Procedure: cppLSip
# What it does:
#   Solves an n by n+1 augmented linear system modulo p using in-place C++ RREF.
# Inputs:
#   - A: integer[8], C_order augmented Matrix; it is destroyed by the solve.
#   - p: prime modulus.
# Outputs:
#   - Returns the unique solution as a list, or FAIL when the matrix is rank deficient or inconsistent.
# Example:
#   x := cppLSip(A,101):
# =============================================================================
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
# =============================================================================
# Procedure: cppLS
# What it does:
#   Non-destructive modular linear-system wrapper for a general augmented Matrix.
# Inputs:
#   - A: n by n+1 augmented Matrix.
#   - p: prime modulus.
# Outputs:
#   - Returns the solution list or FAIL; A is preserved because a machine-integer copy is solved.
# Example:
#   x := cppLS(A,101):
# =============================================================================
cppLS := proc(A::Matrix,p::prime)
    return cppLSip(cppMat64(A,p),p):
end proc:

#  Maple's Modular kernel even below the bound.
MAPLE_FAST_BOUND := 2^31:

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

# =============================================================================
# Procedure: cppEncodeMatrix
# What it does:
#   Encodes a parametric augmented Matrix once as flat sparse monomial data so repeated black-box evaluation can run entirely in C++.
# Inputs:
#   - L: Matrix whose entries are polynomials with rational coefficients in params.
#   - params: ordered list of parameter names.
# Outputs:
#   - Returns table E containing dimensions, sparse term offsets/exponents/raw coefficients, cached reduced coefficients, and reusable work Arrays.
# Example:
#   E := cppEncodeMatrix(<y1*x1+1 | 2>,[y1]):
# =============================================================================
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
# =============================================================================
# Procedure: cppEvalSolve
# What it does:
#   Evaluates an encoded parametric augmented matrix at one parameter point and solves the resulting system modulo p.
# Inputs:
#   - E: encoding returned by cppEncodeMatrix.
#   - point_: list of at least E["nv"] parameter values.
#   - p: prime modulus.
# Outputs:
#   - Returns the solution list for a unique system, or FAIL if the evaluated matrix is singular/inconsistent.
# Example:
#   x := cppEvalSolve(E,[7,11],101):
# =============================================================================
cppEvalSolve := proc(E,point_::list,p::prime)
    local k,rc,nr,nv:

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

# =============================================================================
# Procedure: cppEvalSolveBlock
# What it does:
#   Evaluates and solves an encoded system at a whole block of points in one external call, with output stored transposed by equation.
# Inputs:
#   - E: matrix encoding.
#   - pts: flat point-major Array containing npts parameter points.
#   - npts: number of points.
#   - p: prime modulus.
# Outputs:
#   - Returns a flat Array X whose row i contains x_i across all points, or FAIL if any point is singular/inconsistent.
# Example:
#   X := cppEvalSolveBlock(E,pts,T,101):
# =============================================================================
cppEvalSolveBlock := proc(E,pts::Array,npts::posint,p::prime)
    local k,rc,nr,nv,X:
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
# =============================================================================
# Procedure: cppBlockRow
# What it does:
#   Creates a zero-copy ArrayTools alias for one equation row of the transposed block-solve output.
# Inputs:
#   - X: flat block output Array.
#   - i: 1-based equation index.
#   - npts: physical row stride in X.
#   - m: number of leading entries to expose.
# Outputs:
#   - Returns an Array alias indexed 0..m-1 that shares storage with X.
# Example:
#   Yarr := cppBlockRow(X,2,T,m):
# =============================================================================
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
ROOTS_MODE      := 2:

if USE_POLROOTS64S then
    # =============================================================================
    # External procedure: cppPolRootsext
    # What it does:
    #   Raw external binding to the C++ finite-field polynomial root finder.
    # Inputs:
    #   - d,f: degree and ascending coefficient Array.
    #   - p: prime modulus.
    #   - wsize: scratch-size override.
    #   - outLen,rootsOut: root buffer.
    #   - infoLen,info: metadata buffer.
    # Outputs:
    #   - Returns the C++ status code; info[0] receives the number of roots and rootsOut receives the roots.
    # Example:
    #   rc := cppPolRootsext(d,fA,p,0,d,rts,1,info):
    # =============================================================================
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

    cppPolRootsext := subsop(1=(
                       d,
                       f,
                       p,
                       wsize,
                       outLen,
                       rootsOut,
                       infoLen,
                       info),
                       op(cppPolRootsext)):
    if not type(cppPolRootsext,procedure) then
        error "define_external did not bind cppPolRootsext; check that %1 "
              "exports cppPolRoots",libObj:
    fi:
fi:

#  cppPolRootsOfC(C,p) takes the ASCENDING coefficient list directly -- which
#  is exactly what cppBMM returns -- so no symbolic polynomial is ever built
#  and no coeff() extraction runs.  The old polynomial route cost O(d^2)
#  interpreted work per call: BMEA_poly added d terms into a symbolic sum and
#  cppPolRootsOf then did d+1 coeff() scans of it; at deg 5579 that is ~31M
#  term inspections before the root finder even started.  The C side reduces
#  the coefficients into [0,p) itself.  polroots64s returns DISTINCT roots, so
#  every multiplicity is 1; the shape matches Roots(F) mod p, sorted by root,
#  so R[1][1] = 0 still detects the zero root.
# =============================================================================
# Procedure: cppPolRootsOfC
# What it does:
#   Passes an ascending coefficient list directly to the C++ finite-field root finder, avoiding construction and re-scanning of a symbolic polynomial.
# Inputs:
#   - C: coefficient list [c0,c1,...,cd].
#   - p: prime modulus.
#   - wsize: optional C++ scratch size; 0 requests the default.
# Outputs:
#   - Returns a sorted Maple-style root list [[root,1],...] containing distinct roots modulo p.
# Example:
#   R := cppPolRootsOfC([100,0,1],101):  # roots of x^2-1 mod 101
# =============================================================================
cppPolRootsOfC := proc(C::list,p::prime,wsize::integer := 0)
    local d,fA,rts,info,rc,n,k:
    d := nops(C)-1:
    if d < 1 then
        return []:
    fi:
    fA   := Array(1..d+1,C,datatype=integer[8]):
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

#  Single dispatch point for the three MRFI call sites.  Takes the coefficient
#  list straight from BMEA_coeffs; a symbolic polynomial is built ONLY when the
#  Maple branch actually needs one.
#  n=12 measurements (CECM): 32 bit prime -- Maple Roots 6.00 s vs polroots64s
#  8.90 s through the old polynomial route; 64 bit prime -- polroots64s 21.2 s
#  against a far slower Maple generic path.  Worth re-running ROOTS_MODE=2 at
#  32 bits after this change: the 8.90 included the O(d^2) conversion that no
#  longer exists.
# =============================================================================
# Procedure: rootsMODp
# What it does:
#   Dispatches finite-field root finding between Maple Roots and the C++ polroots64s wrapper according to ROOTS_MODE and prime size.
# Inputs:
#   - C: ascending coefficient list.
#   - Z: variable name used only by the Maple symbolic branch.
#   - p: prime modulus.
# Outputs:
#   - Returns the Maple Roots-style list [[root,multiplicity],...].
# Example:
#   R := rootsMODp([100,0,1],Z,101):
# =============================================================================
rootsMODp := proc(C::list,Z::name,p::prime)
    local i:
    if nops(C) < 2 then
        return []:
    fi:
    if ROOTS_MODE = 1 or not USE_POLROOTS64S then
        return Roots(add(C[i]*Z^(i-1),i=1..nops(C))) mod p:
    fi:
    if ROOTS_MODE = 2 or p >= MAPLE_FAST_BOUND then
        return cppPolRootsOfC(C,p):
    fi:
    return Roots(add(C[i]*Z^(i-1),i=1..nops(C))) mod p:
end proc:

