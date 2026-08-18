kernelopts(numcpus=1):
with(LinearAlgebra):
with(IntegerRelations):

read "./mapleWrapperv2.mpl":
read "./pol_sys.mpl":

ExpectedRawParams := ParametricSystems:-RawParams:
ExpectedParamCount := ParametricSystems:-ParamCount:
DB_system := ParametricSystems:-DBSystem:


get_eqn := proc(Sys,vars)
    option remember:
    return solve(Sys,vars):
end proc:

reording := proc(unordered_soln,num_eqn)
    local component1,ordered_soln,i:
    ordered_soln := [seq(0,i=1..num_eqn)]:
    for i from 1 to num_eqn do
        component1 := get_component(op(1,unordered_soln[i])):
        ordered_soln[component1] := unordered_soln[i]:
    od:
    return ordered_soln:
end proc:

get_component := proc(expression)
    local temp_var:
    temp_var := convert(expression,string):
    return parse(temp_var[2..length(temp_var)]):
end proc:

get_u := proc(M,col,alpha,p)
    local F,U,i;
    F := [seq(convert(M[..,i],list),i=1..col)]:
    U := [seq(cppNewtonInterp(alpha,F[i],x,p),i=1..col)]:
    return U:
end proc:

det_get_u := proc(M,col,alpha,p)
    local F,U:
    F := convert(M[..,col][..numelems(alpha)],list):
    U := cppNewtonInterp(alpha,F,x,p):
    return U:
end proc:

(* Mike's FFGE code. *)

FFGE := proc(A::Matrix, b::Vector, Y::list(name))
local n, m, B, mu, i, j, k, det, x, y, num, r, numterms, f, g, h, mons,
      elim_num_max, elim_post_max, backsub_N_max, y_terms, f_terms, g_terms;
    global ffge_stats;
    numterms := proc(f) if f=0 then 0 elif type(f,`+`) then nops(f) else 1 fi end;
    n,m := op(1,A);
    if n<>m then error "Matrix must be square" fi;
    m := op(1,b);
    if m<>n then error "Matrix and vector must have the same dimension" fi;
    B := <A|b>;
    mu := 1:
    det := 1;
    elim_num_max  := 0:
    elim_post_max := 0:
    for k to n-1 do
        i := k;
        while i<=n and B[i,k]=0 do i := i+1; od;
        if i>n then return 0 fi;
        if i>k then 
            for j from k to n+1 do B[i,j],B[k,j] := B[k,j],B[i,j] od;
            det := -det;
        fi;
        for i from k+1 to n do
            for j from k+1 to n+1 do
                num := expand(B[k,k]*B[i,j]-B[i,k]*B[k,j]);
                divide(num,mu,evaln(B[i,j]));
                if i=k+1 and j=k+1 then
                    #lprint(i,numterms(num),numterms(B[i,j]));
                    if numterms(num)      > elim_num_max  then elim_num_max  := numterms(num)      end if;
                    if numterms(B[i,j])   > elim_post_max then elim_post_max := numterms(B[i,j])   end if;
                fi;
            od;
            B[i,k] := 0;
        od;
        mu := B[k,k];
    od;
    det := det*B[n,n];
    #printf("#det=%d\n",numterms(det));
    # Lipson's back substitution
    y := Vector(n);
    y[n] := B[n,n+1];
    #printf("#y[%d]=%d\n",n,numterms(B[n,n+1]));
    backsub_N_max := 0:
    for i from n-1 by -1 to 1 do
        num := expand(B[i,n+1]*B[n,n]-add(B[i,j]*y[j],j=i+1..n));
        divide(num,B[i,i],evaln(y[i]));
        #printf("#N[%d]=%d  #y[%d]=%d\n",i,numterms(num),i,numterms(y[i]));
        if numterms(num) > backsub_N_max then backsub_N_max := numterms(num) end if;
    od;
    f := Vector(n);
    g := Vector(n);
    for i from 1 to n do
        h := gcd(y[i],B[n,n],evaln(f[i]),evaln(g[i]));
        #printf("#f[%d]=%d #g[%d]=%d #h=%d",i,numterms(f[i]),
             #i,numterms(g[i]),numterms(h));
        #printf("   deg(f[%d])=%d deg(g[%d])=%d\n",i,degree(f[i]),i,degree(g[i]));
    od;
    r := {seq(Y[i]=ithprime(i),i=1..nops(Y))};
    for i from 1 to n do
        coeffs(f[i],Y,'mons');
        mons := subs(r,[mons]);
        m := max(op(mons));
        #printf("f[%d]: max m[%d] = %d = %a\n",i,i,m,ifactor(m));
        coeffs(g[i],Y,'mons');
        mons := subs(r,[mons]);
        m := max(op(mons));
        #printf("g[%d]: max m[%d] = %d = %a\n",i,i,m,ifactor(m));
    od;

    # Stats for the driver.
    y_terms := [seq(numterms(y[i]), i=1..n)]:
    f_terms := [seq(numterms(f[i]), i=1..n)]:
    g_terms := [seq(numterms(g[i]), i=1..n)]:
    ffge_stats := table():
    ffge_stats["elim_num_max"]  := elim_num_max:
    ffge_stats["elim_post_max"] := elim_post_max:
    ffge_stats["backsub_N_max"] := backsub_N_max:
    ffge_stats["y_terms"]       := y_terms:
    ffge_stats["det_terms"]     := numterms(det):
    ffge_stats["f_terms"]       := f_terms:
    ffge_stats["g_terms"]       := g_terms:

    return det,f,g;
end:

Constuct_Sys_Blackbox := proc(Sys,Vars,params)
    local Lin_BB,L,nr,nc,LE:
    L := GenerateMatrix(Sys,Vars,augmented=true):
    print(L);
    nr,nc := op(1,L):

    #  Encode L as a flat sparse monomial table ONCE.  After this the black box
    #  does no symbolic work at all: evaluation at the point and the rref solve
    #  both happen inside a single external call.
    LE := cppEncodeMatrix(L,params):
    Lin_BB := proc(point_::list(integer),p::prime)
        local A,T,subs_values,num_eqn,soln,t0,t_avg,i,j;
        global counter,bb_calls_ndsa,bb_calls_mrfi,t_ndsa_total,t_mrfi_total,bb_phase;
        #uses LinearAlgebra:-Modular:
        counter := counter+1:
        if bb_phase = "NDSA" then bb_calls_ndsa := bb_calls_ndsa+1:
        else                       bb_calls_mrfi := bb_calls_mrfi+1: fi:
        t0 := time():
        #  Only the commented out paths below need subs_values; cppEvalSolve
        #  takes the point as a plain list.
        #subs_values := zip((par,pnt) -> par=pnt,params,point_):
        num_eqn := numelems(Vars):

        #  OLD 32 BIT PATH.  LinearAlgebra:-Modular:-LinearSolve accumulates in
        #  a machine word, so the modulus has to stay below 2^32.
        #
        #A := Mod(p,L,subs_values,integer):
        #T := traperror(LinearSolve(p,A,1)):
        #t_avg := time() - t0:
        #if bb_phase = "NDSA" then t_ndsa_total := t_ndsa_total+t_avg:
        #else                       t_mrfi_total := t_mrfi_total+t_avg: fi:
        #if T = "Matrix is singular." then
        #    return FAIL:
        #fi:
        #soln := convert(A[1..num_eqn,num_eqn+1],list):
        #return soln:

        #  FIRST 64 BIT ATTEMPT.  Correct, but modp(eval(...),p) evaluates over
        #  Z before reducing, so the intermediates are bignums as soon as the
        #  entries of L have any degree in the parameters.  That, not the rref,
        #  is what made this path slow.
        #
        #A := Matrix(nr,nc,(i,j) -> modp(eval(L[i,j],subs_values),p),
        #            datatype=integer[8],order=C_order):
        #soln := cppLSip(A,p):

        #  SECOND 64 BIT ATTEMPT.  cppEvalMatrix64 reduces as it evaluates
        #  instead of afterwards, but Modular:-Mod refuses a 64 bit modulus, so
        #  it falls back to Eval ... mod p, which is still one interpreted
        #  evaluation per matrix entry.  Kept as a fallback, not the hot path.
        #
        #A := cppEvalMatrix64(L,subs_values,p):
        #soln := cppLSip(A,p):

        #  CURRENT 64 BIT PATH.  L was encoded once by cppEncodeMatrix, so the
        #  evaluation at the point AND the rref solve both happen inside one
        #  external call.  Maple only fills the nv parameter values, and no
        #  intermediate ever leaves [0,p).  Good for any prime p < 2^63.
        soln := cppEvalSolve(LE,point_,p):
        t_avg := time() - t0:
        if bb_phase = "NDSA" then t_ndsa_total := t_ndsa_total+t_avg:
        else                       t_mrfi_total := t_mrfi_total+t_avg: fi:
        if soln = FAIL then
            return FAIL:
        fi:
        return soln:
    end proc:
end proc:

#  The whole T x num_var block of points is built in C++.  The arithmetic
#  (beta_[nv-1]*(alpha[np]-s1) + sigma_[nv]) mod p forms products of size up to
#  2^62, which is the edge of Maple's immediate-integer range, so in Maple a
#  large share of these multiplications allocate multi-precision temporaries.
#  In C++ each one is a single 128-bit multiply.  Same name, same signature and
#  same result as before, so the call sites are unchanged.
get_point_on_affine_line := proc(num_var::posint,alpha::list,beta_::list,
                                 sigma_::list,p::prime,T::posint)
    local aArr,bArr,sArr,outArr,i,s:
    global num_lines:
    num_lines := num_lines+1:
    aArr := Array(0..T-1,datatype=integer[8]):
    for i from 1 to T do
        aArr[i-1] := alpha[i]:
    od:
    bArr := Array(0..max(num_var-1,1)-1,datatype=integer[8]):
    for i from 1 to num_var-1 do
        bArr[i-1] := beta_[i]:
    od:
    sArr := Array(0..num_var-1,datatype=integer[8]):
    for i from 1 to num_var do
        sArr[i-1] := sigma_[i]:
    od:
    outArr := Array(0..T*num_var-1,datatype=integer[8]):
    cppAffineLine(T,num_var,aArr,bArr,sArr,p,outArr):
    return [seq([seq(outArr[(s-1)*num_var+i],i=0..num_var-1)],s=1..T)]:
end proc:

MQRFR := proc(r0,r1,t0,t1,p)
    local r,t,q,i,f,g,qmax,lcg;
    r[0] := r0:
    r[1] := r1:
    t[0] := t0: 
    t[1] := t1:
    f := r0:
    g := t1: 
    qmax := 0: 
    i := 1:
    while r[i] <> 0 do
        q[i] := Quo(r[i-1],r[i],x,'r[i+1]') mod p:
        if degree(q[i],x)>qmax then
            qmax := degree(q[i],x):
            f := r[i]:
            g := t[i]:
        fi:
        t[i+1] := Expand(t[i-1]-q[i]*t[i]) mod p:
        i := i+1:
        if qmax <= 1 or gcd(f,g) <> 1 or g = 0 then 
            FAIL: 
        fi:
    od:
    lcg := lcoeff(g):
    return f/lcg mod p,g/lcg mod p,qmax,lcg:
end proc:

# HFTRFR implementation

HFTRFR := proc(r0, r1, t0, t1, p)
    local r, t, q, i, f, g, qmax, Lambda, fc, gc, ilc, rts, badset:
    r[0] := r0:  r[1] := r1:
    t[0] := t0:  t[1] := t1:
    f := r0:  g := t1:  qmax := 0:  i := 1:
    while r[i] <> 0 do
        q[i] := Quo(r[i-1], r[i], x, 'r[i+1]') mod p:
        if degree(q[i], x) > qmax then
            qmax := degree(q[i], x):
            f := r[i]:
            g := t[i]:
        fi:
        t[i+1] := Expand(t[i-1] - q[i]*t[i]) mod p:
        i := i+1:
    od:
    Lambda := Gcd(f, g) mod p:              # bad polynomial
    fc := Quo(f, Lambda, x) mod p:
    gc := Quo(g, Lambda, x) mod p:
    ilc := 1/lcoeff(gc, x) mod p:           # make gc monic
    fc := fc*ilc mod p:
    gc := gc*ilc mod p:
    if degree(Lambda, x) = 0 then
        badset := []:
    else
        rts := Roots(Lambda) mod p:
        badset := sort([seq(rts[i][1], i=1..nops(rts))]):
    fi:
    return fc, gc, degree(fc, x), degree(gc, x), qmax, badset, Lambda:
end proc:

# DFTRFR implementation

DFTRFR := proc(M, U, degN, degD, E, p)
    local Rp, R, Tp, T, Q, d, num, den, ilc, cg:
    if degree(M, x) <= degN + degD + 2*E then
        return FAIL:                                   
    fi:
    if U = 0 then
        return 0:                                      
    fi:
    cg := Gcd(M, U) mod p:

    # GCD rejection test

    if degree(cg, x) > degN + E then
        return FAIL:                                   
    fi:
    Rp := M:  R := U:
    Tp := 0:  T := 1:
    while R <> 0 and degree(T, x) <= degD + E do
        Q      := Quo(Rp, R, x) mod p:
        Rp, R  := R, (Expand(Rp - Q*R) mod p):
        Tp, T  := T, (Expand(Tp - Q*T) mod p):
    od:
    d   := Gcd(Rp, Tp) mod p:
    num := Quo(Rp, d, x) mod p:
    den := Quo(Tp, d, x) mod p:
    if degree(num, x) > degN or degree(den, x) > degD then
        return FAIL:                                   
    fi:
    ilc := 1/lcoeff(den, x) mod p:
    return (num*ilc mod p) / (den*ilc mod p):
end proc:

# Corrupt values

inject_faults := proc(Y::list, E::nonnegint, p::prime)
    local n, chosen, Yc, i, delta:
    n := numelems(Y):
    if E = 0 or n = 0 then return Y, []: fi:
    if E >= n then error "cannot inject %1 faults into %2 points", E, n: fi:
    chosen := combinat:-randcomb(n, E):     
    Yc := Y:
    for i in chosen do
        delta := rand(1..p-1)():
        Yc := subsop(i = ((Y[i] + delta) mod p), Yc):
    od:
    return Yc, sort([op(chosen)]):
end proc:

get_u_faulty := proc(M,col,alpha,p)
    local F,U,i,Fi;
    global FAULT_E;
    F := [seq(convert(M[..,i],list),i=1..col)]:
    U := []:
    for i from 1 to col do
        Fi := inject_faults(F[i][1..numelems(alpha)],FAULT_E,p)[1]:
        U := [op(U),cppNewtonInterp(alpha,Fi,x,p)]:
    od:
    return U:
end proc:

#  Berlekamp-Massey runs in C++ (cppBMM -> cppBM).  It is split in two here:
#  BMEA_coeffs returns the raw coefficient list, and BMEA_poly turns that list
#  into a polynomial in Z.  Only the degree is needed to test convergence, and
#  the degree is just nops(L)-1, so the (large) symbolic polynomial is built
#  only on the rounds where Roots is actually called.
BMEA_coeffs := proc(v::{Vector,list}, p::posint)
    if numelems(v)=0 then
        return []:
    fi:
    return cppBMM(v,p):
end proc:

BMEA_poly := proc(L::list, Z::name)
    local d,i:
    if L = [] then
        return 1:
    fi:
    d := nops(L)-1:
    return add(L[i+1]*Z^i,i=0..d):
end proc:

#  Kept for callers that want the polynomial directly.
BMEA := proc(v::{Vector,list}, p::posint, Z::name)
    return BMEA_poly(BMEA_coeffs(v,p),Z):
end proc:

generate_monomials := proc(roots_,num_var,prime_points,vars)
    local m,mm,i,j,counter_,M_,rem:
    M_ := Vector(numelems(roots_),0):
    for i from 1 to numelems(roots_) do
        if roots_[i] = 0 then 
            return FAIL: 
        fi:
        mm := roots_[i]:
        m := 1:
        for j from 1 to numelems(prime_points) do
            counter_ := 0:
            while mm mod prime_points[j] = 0 do
                mm := iquo(mm,prime_points[j],'rem'):
                counter_ := counter_+1:
            od:
            m := m*vars[j]^counter_:
        od:
        M_[i] := m:
        if mm <> 1 then
            print("Warning: residue mm=", mm, " (should be 1) for root[", i, "]"):
            return FAIL:
        fi:
    od:
    return convert(M_,list):
end proc:

Zippel_Transpose_Vandermonde_solver := proc(y::{Vector,list}, terms::integer,
                                            roots_::list, lambda_::polynom,
                                            p::integer)
    if terms = 0 then 
    return []: 
    fi:
    return cppVS(y[1..terms],roots_[1..terms],p,1):
end proc:

construct_final_polynomial := proc(coeff_,Monomials)
    local i,f:
    f := add(coeff_[i]*Monomials[i],i=1..numelems(coeff_)):
end proc:

NDSA := proc(B,sigma_,beta_,num_var,p,num_points,num_eqn)
    local correct_degree,T,alpha,m,Psi_alpha,Y,u,dq,i,r,
          lin_sys,temp,result,count,M,row,col,DQ,MQRFR_done,
          t_helper:
    global t_cppNewton_total,t_mqrfr_total,t_NDSA_comp:
    print("In NDSA"):
    MQRFR_done := [seq(false,i=1..num_eqn)]:
    correct_degree := false:
    lin_sys := false:
    T := num_points:
    temp := []:
    result := [seq([],i=1..num_eqn)]:
    count := 0:
    while not correct_degree do
        count := count+1:
        print("NDSA: T = ", T):
        r := rand(p):
        t_helper := time():
        alpha := [seq(r(),i=1..T)]:
        m := Expand(product(x-alpha[j],j=1..T)) mod p:
        Psi_alpha := get_point_on_affine_line(num_var,alpha,beta_,sigma_,p,T):
        print(Psi_alpha):
        t_NDSA_comp := t_NDSA_comp+(time()-t_helper):
        Y := [seq(B(Psi_alpha[i],p),i=1..T)]:
        M := Matrix(Y):
        print(Y):
        print(M):
        row,col := Dimension(M):
        if row = 1 then
            lin_sys := false:
            if FAULT_ON then Y := inject_faults(Y,FAULT_E,p)[1]: fi:
            t_helper := time():
            u := cppNewtonInterp(alpha,Y,x,p):
            t_cppNewton_total := t_cppNewton_total+(time()-t_helper):
            t_helper := time():
            result := [[cppHFTRFR(m,u,x,p)]]:
            t_mqrfr_total := t_mqrfr_total +(time()-t_helper):
            dq := result[1][5]:                 
        else
            lin_sys := true:
            t_helper := time():
            if FAULT_ON then
                u := get_u_faulty(M,col,alpha,p):
            else
                u := get_u(M,col,alpha,p):
            fi;
            t_cppNewton_total := t_cppNewton_total+(time()-t_helper):
        fi:
        if lin_sys then
            t_helper := time():
            for i from 1 to nops(u) do
                if MQRFR_done[i] then 
                    next:
                fi:
                temp := [op(temp),cppHFTRFR(m,u[i],x,p)]:
                result[i] := temp:
                temp := []:
            od:
            t_mqrfr_total := t_mqrfr_total+(time()-t_helper):
            DQ := [seq(result[i][5],i=1..nops(result))];   
            for i from 1 to numelems(DQ) do
                if DQ[i] > 1 then 
                    MQRFR_done[i] := true:
                fi:
            od:
            dq := min(op(DQ)):
        fi:
        if dq > 1 then
            print("NDSA: Termination condition met."):
            return result,lin_sys:
        else
            print("NDSA: MQRFR failed -> Doubling T value"):
            T := T*2:
            for i from 1 to num_eqn do
                if not MQRFR_done[i] then
                    result[i] := []:
                fi:
            od:
            DQ := []:
        fi:
        if T > 2^20 then
            print("NDSA: Hit safety break at T = ", T):
            return result,lin_sys:
        fi:
    od:
end proc:

MRFI := proc(B,num_vars::integer,num_eqn::integer,vars::list,p::integer)
    local i,j,k,s,Primes,direction,sigma_,num_eval,den_eval,u,mon,
          num_eval_n,den_eval_n,
          numerator_done,denominator_done,Tcur,jDone,
          mqrfr_results,lin_sys,num_points_mqrfr,
          Numerators,Denominiators,deg_num,deg_den,
          lambda_num,lambda_den,terms_num,terms_den,
          R_num,R_den,Roots_num_eval,Roots_den_eval,
          num_mono,den_mono,coeff_num,coeff_den,
          final_num,final_den,temp,common_den_flag,
          bmea_done,temp_den,all_den_done,all_done,max_num_points_mqrfr,
          init_sigma,sampleCounts,mMax,tempG,
          sigma_j,sigma_run,alphaVal,
          Psi_alpha,BBvals,m_i,
          alphaArr,Yarr,outArr,Yfault,
          evalCap,newCap,tmpV,lam_c,
          gc0,alloc0,t_helper,tempR:
    global mrfi_stats,bb_phase,
           t_cppNewton_total,t_cppRR_total,t_bmea_total,t_zippel_total,
           t_roots_total,t_genmono_total,t_construct_total,
           t_points_total,t_seqeval_total,t_normalize_total:

    gc0    := kernelopts(gctimes):
    alloc0 := kernelopts(bytesused):

    Primes := [seq(ithprime(i),i=1..num_vars)]:
    tempR := rand(0..p-1):
    direction := [seq(tempR(),i=1..num_vars)]:
    init_sigma := Primes:
    sigma_ := []:
    num_eval := table(): 
    num_eval_n := table():
    den_eval := table(): 
    den_eval_n := table():
    numerator_done   := [seq(false,i=1..num_eqn)]:
    denominator_done := [seq(false,i=1..num_eqn)]:
    bmea_done := [seq(false,i=1..num_eqn)]:
    all_done := true:

    lambda_num := table(): 
    terms_num := table(): 
    R_num := table():
    lambda_den := table(): 
    terms_den := table(): 
    R_den := table():
    final_num  := table(): 
    final_den := table():
    num_mono   := table(): 
    coeff_num := table():
    den_mono   := table(): 
    coeff_den := table():

    for i from 1 to num_eqn do
        lambda_num[i] := []: 
        terms_num[i] := []: 
        R_num[i] := []:
        lambda_den[i] := []: 
        terms_den[i] := []: 
        R_den[i] := []:
    od:

    Tcur := 4:
    num_points_mqrfr := [seq(0,i=1..num_eqn)]:

    bb_phase := "NDSA":
    mqrfr_results,lin_sys := NDSA(B,init_sigma,direction,
                                    num_vars,p,Tcur,num_eqn):
    print(mqrfr_results):
    print(lin_sys):
    bb_phase := "MRFI":

    Numerators := [seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:
    Denominiators := [seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
    
    # Optimization from previous implementations.
    # Evaluations live in machine-integer Vectors, not tables: an append is a
    # word store instead of a hash insert, and BMEA/Vandermonde can be handed
    # a contiguous slice instead of a freshly built list.

    evalCap := 16:
    for k from 1 to num_eqn do
        num_eval[k] := Vector[column](evalCap,datatype=integer[8]):
        num_eval_n[k] := 0:
        den_eval[k] := Vector[column](evalCap,datatype=integer[8]):
        den_eval_n[k] := 0:
    od:

    deg_num := [seq(degree(Numerators[i],x),i=1..nops(Numerators))]:
    deg_den := [seq(degree(Denominiators[i],x),i=1..nops(Denominiators))]:
    for i from 1 to numelems(deg_den) do
        num_points_mqrfr[i] := deg_num[i]+deg_den[i]+1+2*FAULT_E:   # +2E for DFTRFR
    od:
    max_num_points_mqrfr := max(op(deg_num))+max(op(deg_den))+1+2*FAULT_E:

    sampleCounts := num_points_mqrfr:
    mMax := max_num_points_mqrfr:
    tempG := rand(1..p-1):
    jDone := 0:
    sigma_run := [seq(Primes[k]^jDone mod p,k=1..nops(Primes))]:

    common_den_flag := true:
    for k from 2 to num_eqn do
        if Denominiators[k] <> Denominiators[1] then
            common_den_flag := false:
            break:
        fi:
    od:

    alphaVal := [seq(tempG(),s=1..mMax)]:
    alphaArr := Array(0..mMax-1,datatype=integer[8]):
    for s from 1 to mMax do
        alphaArr[s-1] := alphaVal[s]:
    od:
    Yarr   := Array(0..mMax-1,datatype=integer[8]):
    outArr := Array(0..1,datatype=integer[8]):

    while true do
        if evalCap < 2*Tcur then
            newCap := 2*Tcur:
            for k from 1 to num_eqn do
                tmpV := Vector[column](newCap,datatype=integer[8]):
                if num_eval_n[k] > 0 then
                    tmpV[1..num_eval_n[k]] := num_eval[k][1..num_eval_n[k]]:
                fi:
                num_eval[k] := tmpV:
                tmpV := Vector[column](newCap,datatype=integer[8]):
                if den_eval_n[k] > 0 then
                    tmpV[1..den_eval_n[k]] := den_eval[k][1..den_eval_n[k]]:
                fi:
                den_eval[k] := tmpV:
            od:
            evalCap := newCap:
        fi:
        for j from jDone+1 to 2*Tcur do
            sigma_run := [seq(sigma_run[k]*Primes[k] mod p,k=1..nops(Primes))]:
            sigma_j := sigma_run:
            if lin_sys then
                t_helper := time():
                Psi_alpha := get_point_on_affine_line(num_vars,alphaVal,
                                                     direction,sigma_j,
                                                     p,mMax):
                t_points_total := t_points_total+(time()-t_helper):
                BBvals := [seq(B(Psi_alpha[s], p),s=1...mMax)]:
                for i from 1 to num_eqn do

                    # Optimization 2:
                    # Skip work that cannot be used.  When the denominator is
                    # common, only equation 1 supplies the denominator
                    # sequence (the assembly copies coeff_den[1] to the rest),
                    # so once equation i>1 has its numerator it has nothing
                    # left to contribute and is skipped entirely.
                    if numerator_done[i] and
                       (denominator_done[i] or (common_den_flag and i > 1)) then
                        next:
                    fi:
                    m_i := sampleCounts[i]:
                    for s from 1 to m_i do
                        Yarr[s-1] := BBvals[s][i]:
                    od:
                    if FAULT_ON then
                        Yfault := inject_faults([seq(Yarr[s-1],s=1..m_i)],
                                                FAULT_E,p)[1]:
                        for s from 1 to m_i do
                            Yarr[s-1] := Yfault[s]:
                        od:
                    fi:

                    #  All interpolation + build modulus + DFTRFR + evaluate at
                    #  sigma are all done inside C++ returning two machine integers.

                    t_helper := time():
                    if cppFTREval(m_i,alphaArr,Yarr,sigma_j[1],
                                  deg_num[i],deg_den[i],FAULT_E,p,
                                  outArr) = FAIL then
                        error "cppFTREval failed at eqn %1 "
                              "(increase points or lower FAULT_E)", i:
                    fi:
                    t_cppRR_total := t_cppRR_total+(time()-t_helper):
                    if not numerator_done[i] then
                        num_eval_n[i] := num_eval_n[i]+1:
                        num_eval[i][num_eval_n[i]] := outArr[0]:
                    fi:
                    if (not denominator_done[i]) and
                       (i = 1 or not common_den_flag) then
                        den_eval_n[i] := den_eval_n[i]+1:
                        den_eval[i][den_eval_n[i]] := outArr[1]:
                    fi:
                od:
            else
                bb_phase := "NDSA":
                mqrfr_results,lin_sys := NDSA(B,sigma_j,direction,
                                                num_vars,p,
                                                max_num_points_mqrfr,1):
                bb_phase := "MRFI":
                Numerators := [seq(mqrfr_results[k][1],
                                   k=1..nops(mqrfr_results))]:
                Denominiators := [seq(mqrfr_results[k][2],
                                      k=1..nops(mqrfr_results))]:
                for k from 1 to num_eqn do
                    t_helper := time():
                    if not numerator_done[k] then
                        num_eval_n[k] := num_eval_n[k]+1:
                        num_eval[k][num_eval_n[k]] := Eval(Numerators[k],x=sigma_j[1]) mod p:
                    fi:
                    if not denominator_done[k] then
                        den_eval_n[k] := den_eval_n[k]+1:
                        den_eval[k][den_eval_n[k]] := Eval(Denominiators[k],x=sigma_j[1]) mod p:
                    fi:
                    t_seqeval_total := t_seqeval_total+(time()-t_helper):
                od:
            fi:
        od:

        jDone := 2*Tcur:

        all_done := true:
        for k from 1 to num_eqn do
            if numerator_done[k] then 
                next:
            fi:
            t_helper := time():
            lam_c := BMEA_coeffs(num_eval[k][1..num_eval_n[k]],p):
            t_bmea_total := t_bmea_total+(time()-t_helper):
            terms_num[k] := `if`(lam_c = [], 0, nops(lam_c)-1):
            if terms_num[k] < iquo(num_eval_n[k], 2) then
                lambda_num[k] := BMEA_poly(lam_c,Z):
                t_helper := time():
                #R_num[k] := Roots(lambda_num[k]) mod p:
                R_num[k] := cppRootsOf(lambda_num[k],Z,p):
                t_roots_total := t_roots_total+(time()-t_helper):
                if R_num[k] <> [] then
                    if nops(R_num[k]) > 0 and R_num[k][1][1] = 0 then
                        R_num[k] := remove(x -> x=[0,1],R_num[k]):
                        terms_num[k] := terms_num[k]-1:
                    fi:
                    if nops(R_num[k]) = terms_num[k] then
                        numerator_done[k] := true:
                    fi:
                fi:
            fi:
        od:

        all_den_done := true:
        if common_den_flag then
            if not denominator_done[1] then
                t_helper := time():
                lam_c := BMEA_coeffs(den_eval[1][1..den_eval_n[1]],p):
                t_bmea_total := t_bmea_total+(time()-t_helper):
                terms_den[1] := `if`(lam_c = [], 0, nops(lam_c)-1):
                if terms_den[1] < iquo(den_eval_n[1],2) then
                    lambda_den[1] := BMEA_poly(lam_c,Z):
                    t_helper := time():
                    #R_den[1] := Roots(lambda_den[1]) mod p:
                    R_den[1] := cppRootsOf(lambda_den[1],Z,p):
                    t_roots_total := t_roots_total+(time()-t_helper):
                    if R_den[1] <> [] then
                        if nops(R_den[1]) > 0 and R_den[1][1][1] = 0 then
                            R_den[1] := remove(x -> x=[0,1],R_den[1]):
                            terms_den[1] := terms_den[1]-1:
                        fi:
                        if nops(R_den[1]) = terms_den[1] then
                            for k from 1 to num_eqn do
                                denominator_done[k] := true:
                            od:
                        fi:
                    fi:
                fi:
            fi:
        else
            for k from 1 to num_eqn do
                if denominator_done[k] then next; end if;
                t_helper := time():
                lam_c := BMEA_coeffs(den_eval[k][1..den_eval_n[k]],p):
                t_bmea_total := t_bmea_total+(time()-t_helper):
                terms_den[k] := `if`(lam_c = [], 0, nops(lam_c)-1):
                if terms_den[k] < iquo(den_eval_n[k],2) then
                    lambda_den[k] := BMEA_poly(lam_c,Z):
                    t_helper := time():
                    #R_den[k] := Roots(lambda_den[k]) mod p:
                    R_den[k] := cppRootsOf(lambda_den[k],Z,p):
                    t_roots_total := t_roots_total+(time()-t_helper):
                    if R_den[k] <> [] then
                        if nops(R_den[k]) > 0 and R_den[k][1][1] = 0 then
                            R_den[k] := remove(x -> x=[0,1],R_den[k]):
                            terms_den[k] := terms_den[k]-1:
                        fi:
                        if nops(R_den[k]) = terms_den[k] then
                            denominator_done[k] := true:
                        fi:
                    fi:
                fi:
            od:
        fi:

        for i from 1 to num_eqn do
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        od:
        all_done := true;
        for i from 1 to num_eqn do 
            all_done := all_done and bmea_done[i]:
        od:
        if all_done then 
            break: 
        fi:

        Tcur := 2*Tcur:
        if Tcur > 2^20 then
            print("MRFI: Hit safety break at Tcur=",Tcur):
            break:
        fi:
    od:

    Roots_num_eval := [seq([seq(r[1],r in R_num[k])],k=1..num_eqn)]:
    if common_den_flag then
        Roots_den_eval := [[seq(r[1],r in R_den[1])]]:
        for k from 2 to num_eqn do
            Roots_den_eval := [op(Roots_den_eval),Roots_den_eval[1]]:
        od:
    else
        Roots_den_eval := [seq([seq(r[1],r in R_den[k])],k=1..num_eqn)]:
    fi:

    for k from 1 to num_eqn do
        t_helper := time():
        temp := generate_monomials(Roots_num_eval[k],num_vars,Primes,vars):
        t_genmono_total := t_genmono_total+(time()-t_helper):
        if temp = FAIL then 
            return FAIL: 
        fi:
        num_mono[k] := temp:
        t_helper := time():
        coeff_num[k] := Zippel_Transpose_Vandermonde_solver(num_eval[k][1..terms_num[k]],
                              terms_num[k],Roots_num_eval[k],lambda_num[k],p):
        t_zippel_total := t_zippel_total+(time()-t_helper):
    od:

    if common_den_flag then
        t_helper := time():
        temp := generate_monomials(Roots_den_eval[1],num_vars,Primes,vars):
        t_genmono_total := t_genmono_total+(time()-t_helper):
        if temp = FAIL then 
            return FAIL:
        fi:
        den_mono[1] := temp:
        t_helper := time():
        coeff_den[1] := Zippel_Transpose_Vandermonde_solver(den_eval[1][1..terms_den[1]],
                              terms_den[1],Roots_den_eval[1],lambda_den[1],p):
        t_zippel_total := t_zippel_total+(time()-t_helper):
        for k from 2 to num_eqn do
            den_mono[k] := den_mono[1]:
            coeff_den[k] := coeff_den[1]:
        od:
    else
        for k from 1 to num_eqn do
            t_helper := time():
            temp := generate_monomials(Roots_den_eval[k],num_vars,Primes,vars):
            t_genmono_total := t_genmono_total+(time()-t_helper):
            if temp = FAIL then 
                return FAIL:
            fi:
            den_mono[k] := temp:
            t_helper := time():
            coeff_den[k] := Zippel_Transpose_Vandermonde_solver(den_eval[k][1..terms_den[k]],
                              terms_den[k],Roots_den_eval[k],lambda_den[k],p):
            t_zippel_total := t_zippel_total+(time()-t_helper):
        od:
    fi:

    for k from 1 to num_eqn do
        t_helper := time():
        lcoeff(add(mon, mon in den_mono[k]),vars,'mon'):
        if not member(mon, den_mono[k], 'i') then
            error "Bug in leading monomial":
        fi:
        u := 1/coeff_den[k][i] mod p:
        coeff_num[k] := u*coeff_num[k] mod p:
        coeff_den[k] := u*coeff_den[k] mod p:
        t_normalize_total := t_normalize_total+(time()-t_helper):
        t_helper := time():
        final_num[k] := construct_final_polynomial(coeff_num[k],num_mono[k]):
        final_den[k] := construct_final_polynomial(coeff_den[k],den_mono[k]):
        t_construct_total := t_construct_total+(time()-t_helper):
    od:

    mrfi_stats := table():
    mrfi_stats["terms_num"] := [seq(terms_num[k],k=1..num_eqn)]:
    mrfi_stats["terms_den"] := [seq(terms_den[k],k=1..num_eqn)]:
    mrfi_stats["deg_num"] := deg_num:
    mrfi_stats["deg_den"] := deg_den:
    mrfi_stats["common_den_flag"] := common_den_flag:
    mrfi_stats["mMax"] := mMax:
    mrfi_stats["sampleCounts"] := sampleCounts:
    mrfi_stats["gc_count"] := kernelopts(gctimes)-gc0:
    mrfi_stats["gc_Mbytes"] := (kernelopts(bytesused)-alloc0)/1048576.0:

    return final_num,final_den:
end proc:

RandRational := proc(N::posint)
    return proc() 
    local a,b:
        a := rand(-N..N)():
        b := rand(1..N)():
        if a = 0 then 
            0: 
        else 
            a/b: 
        fi:
    end proc:
end proc:

get_data := proc(test_case)
    local Sys,Vars,i,ii,jj,params,ff,gg,n:
    if nargs=1 then
        if test_case = "bspline" then
            Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1,
                    x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17,
                    -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
                    -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2,
                    -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13,
                    -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
                    2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21,
                    2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7,
                    2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9,
                    -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
                    2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16,
                    2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
        elif test_case = "small_sys_low_deg" then
            Sys := {x1+y1*x2+y2*x3-1, y2*x1+x2+y1*x3-2, (y1-y2)*x1-x2+y2*x3-7}:
        elif test_case = "small_Sys" then
            Sys := {x1+y1*x2+y1-3, y2*x1+x2+y1-1}:
        elif test_case = "mike" then
            Sys := {y1*x1+y1*x2-1, y1*y2*x1-x2-1}:
        elif test_case = "example" then
            Sys := {(y1*y2-1)*x1 + (y1^2-2*y1+3),
                    (y1*y2-1)*x2 + (y1*y2-y1-3*y2+1)}:
        elif test_case = "bsbug" then
            Sys := {(2*y3^2*y4 - y3*y4^2 + 3*y3*y4 - y4^2 + y3 + y4 + 1)*x1
                    = y3*y4^2}:
        elif test_case = 1 then
            ff := y:
            gg := x-4:
            Vars := [x, y]:
            return Vars,ff,gg,numelems(Vars),1,Vars:
        fi:
    elif nargs>1 then
        Vars := [seq(x||i,i=1..args[2])]:
        if test_case = "rand" then
            ff := randpoly(Vars,terms=args[3]):
            gg := randpoly(Vars,terms=args[4]):
            return Vars,ff,gg,numelems(Vars),1,Vars:
        elif test_case = "rat_rand" then
            ff := randpoly(Vars,coeffs=RandRational(args[5]),terms=args[3]):
            gg := randpoly(Vars,coeffs=RandRational(args[6]),terms=args[4]):
            return Vars,ff,gg,numelems(Vars),1,Vars:
        elif test_case = "TP" then
            n := args[2]:
            Vars := [seq(x||i,i=1..n)]:
            params := [seq(y||i,i=1..n)]:
            Sys := [seq(add(y||(abs(ii-jj)+1)*x||jj,jj=1..n)-1,
                            ii=1..n)]:
            return Sys,Vars,params,numelems(params),numelems(Vars):
        fi:
    fi:
    Vars := [seq(x||i,i=1..nops(Sys))]:
    params := convert(indets(Sys) minus convert(Vars,set),list):
    return Sys,Vars,params,nops(params),nops(Vars):
end proc:

# ---- system selection ------------------------------------------------------
# Change SYSTEM_ID to any ID printed by ParametricSystems:-List().
# All matrix construction, RHS construction, physical parameter lists, and
# y1,y2,... remapping are fetched from parametric_systems.mpl.
#
# Examples: "S1", "R2", "P1", "P13", "P35", "P40".
SYSTEM_ID := "S1":
if not type(SYSTEM_ID, string) then SYSTEM_ID := convert(SYSTEM_ID, string): fi:

# Freeze the selected family for this benchmark run.
RUN_SYSTEM_ID := SYSTEM_ID:

#test_prime := prevprime(2^63-1):
test_prime := prevprime(2^31-1):

# n is the scalable input knob used by the selected family.
# For q-by-q grid systems q=n; for P40 n is the number of QBD levels.
n_min := 4:
n_max := 8:
do_verify := false:
do_ffge := false:
summary := []:

# Controls for FTRFR:

# FAULT_ON = false  reproduces the fault free pipeline (HFTRFR/DFTRFR with
#                   E=0 behave like MQRFR/ordinary RR).
# FAULT_ON = true   injects FAULT_E faults per univariate reconstruction.
#                   HFTRFR discovers degrees + bad set and DFTRFR reconstructs.

FAULT_E  := 0:
FAULT_ON := false:

#  Garbage collector optimization:

#  MRFI allocates heavily in short bursts.  The default collection threshold
#  makes Maple sweep very frequently at large n; raising it trades memory for
#  far fewer collections.  Lower this if the machine runs short of RAM.
kernelopts(gcfreq = 32*10^6):

for n_test from n_min to n_max do
    Sys,Vars,params,num_vars,num_eqn := DB_system(RUN_SYSTEM_ID,n_test):

    

    # Validate every catalog system before any black-box calls are made.
    expected_num_vars := ExpectedParamCount(RUN_SYSTEM_ID,n_test):
    if num_vars <> expected_num_vars then
        error "%1 parameter mismatch at n=%2: expected %3, got %4",
              RUN_SYSTEM_ID, n_test, expected_num_vars, num_vars:
    fi:
    printf("Parameter validation: PASS  (expected=%d, actual=%d)\n",
           expected_num_vars, num_vars):

    printf("\n==== SYSTEM %s : %s   (n = %d, params = %d) ====\n",
           RUN_SYSTEM_ID,
           ParametricSystems:-Name(RUN_SYSTEM_ID),
           n_test, num_vars):

    print(Sys);
    print(Vars);
    print(params);
    print(num_vars);
    print(num_eqn);

    counter := 0:
    num_lines := 0:
    bb_phase := "MRFI":
    bb_calls_ndsa := 0:
    bb_calls_mrfi := 0:
    t_ndsa_total := 0.0:       
    t_mrfi_total := 0.0:
    t_cppNewton_total := 0.0:  
    t_mqrfr_total := 0.0:
    t_cppRR_total := 0.0:      
    t_bmea_total := 0.0:
    t_zippel_total := 0.0:     
    t_roots_total := 0.0:
    t_genmono_total := 0.0:    
    t_construct_total := 0.0:
    t_points_total := 0.0:     
    t_seqeval_total := 0.0:
    t_normalize_total := 0.0:  
    t_iratrecon_total := 0.0:
    t_NDSA_comp := 0.0:        
    t_mrfi_wall := 0.0:
    t_ffge := 0.0:

    B := Constuct_Sys_Blackbox(Sys,Vars,params):
    Num,Den := 'Num','Den':
    status := "OK":
    mrfi_out := []:
    t_mrfi_wall_start := time():
    try
        mrfi_out := [MRFI(B,num_vars,num_eqn,params,test_prime)]:
    catch:
        status := cat("ERROR: ",StringTools:-FormatMessage(lastexception[2..-1])):
        printf("MRFI threw: %s\n",status):
    end try:
    t_mrfi_wall := time()-t_mrfi_wall_start:
    mrfi_calls := counter:
    if status = "OK" then
        if nops(mrfi_out) = 2 and mrfi_out <> [FAIL, FAIL] then
            Num := mrfi_out[1]:
            Den := mrfi_out[2]:
        else
            status := "MRFI returned FAIL":
        fi:
    fi:

    if status = "OK" then
        stats_terms_num := mrfi_stats["terms_num"]:
        stats_terms_den := mrfi_stats["terms_den"]:
        stats_deg_num := mrfi_stats["deg_num"]:
        stats_deg_den := mrfi_stats["deg_den"]:
        stats_common := mrfi_stats["common_den_flag"]:
        stats_mMax := mrfi_stats["mMax"]:
        stats_gc_count := mrfi_stats["gc_count"]:
        stats_gc_Mbytes := mrfi_stats["gc_Mbytes"]:
    else
        stats_terms_num := []: 
        stats_terms_den := []:
        stats_deg_num := []: 
        stats_deg_den := []:
        stats_common := false: 
        stats_mMax := 0:
        stats_gc_count := 0: 
        stats_gc_Mbytes := 0.0:
    fi:

    ffge_run := do_ffge:
    if ffge_run then
        A_ffge,b_ffge := GenerateMatrix(Sys,Vars):
        Y_ffge := params:
        ffge_status := "OK":
        try
            det_ffge,f_ffge,g_ffge := FFGE(A_ffge,b_ffge,Y_ffge):
        catch:
            ffge_status := cat("ERROR: ",
                               StringTools:-FormatMessage(lastexception[2..-1])):
            print("FFGE threw: ",ffge_status):
        end try:

        if ffge_status = "OK" then
            ffge_terms_num := ffge_stats["f_terms"]:
            ffge_terms_den := ffge_stats["g_terms"]:
            ffge_y_terms := ffge_stats["y_terms"]:
            ffge_det_terms := ffge_stats["det_terms"]:
            ffge_elim_swell := ffge_stats["elim_num_max"]:
            ffge_backsub_swell := ffge_stats["backsub_N_max"]:
            # ffge_bench_calls := 1:
            # t_ffge_start := time():
            # to ffge_bench_calls do
            #     FFGE(A_ffge,b_ffge,Y_ffge):
            # od:
            # t_ffge := evalf((time()-t_ffge_start)/ffge_bench_calls):
        else
            ffge_terms_num := []:  
            ffge_terms_den := []:
            ffge_y_terms := []:  
            ffge_det_terms := -1:
            ffge_elim_swell := -1:  
            ffge_backsub_swell := -1:
            t_ffge := 0.0:
        fi:
    else
        ffge_status := `if`(do_ffge, "SKIPPED (rational entries)", "SKIPPED"):
        t_ffge := 0.0:
        ffge_terms_num := []:  
        ffge_terms_den := []:
        ffge_y_terms := []:  
        ffge_det_terms := -1:
        ffge_elim_swell := -1:  
        ffge_backsub_swell := -1:
    fi:

    if status = "OK" then
        Ratrecon_num := table():
        Ratrecon_den := table():
        Final_rat_poly := table():
        t_helper := time():
        for i from 1 to num_eqn do
            Ratrecon_num[i] := iratrecon(Num[i],test_prime):
            Ratrecon_den[i] := iratrecon(Den[i],test_prime):
            Final_rat_poly[i] := Ratrecon_num[i]/Ratrecon_den[i]:
        od:
        t_iratrecon_total := time()-t_helper:
        solution_list := [seq(Final_rat_poly[i], i=1..num_eqn)]:

        if do_verify then
            try
                A_ref,b_ref := GenerateMatrix(Sys,Vars):
                x_ref := LinearSolve(A_ref,b_ref):
                all_match := true:
                for i from 1 to num_eqn do
                    diff_i := normal(Final_rat_poly[i]-x_ref[i]):
                    if diff_i <> 0 then
                        all_match := false:
                    fi:
                od:
                if all_match then
                    printf("PASS  n=%2d\n", n_test):
                    summary := [op(summary),
                                [n_test, "PASS",mrfi_calls,
                                 bb_calls_ndsa,bb_calls_mrfi,
                                 stats_terms_num,stats_terms_den,
                                 stats_deg_num,stats_deg_den,
                                 stats_common,stats_mMax,
                                 ffge_status,t_ffge,
                                 ffge_terms_num,ffge_terms_den,
                                 ffge_y_terms,ffge_det_terms,
                                 ffge_elim_swell,ffge_backsub_swell,
                                 t_ndsa_total,t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total,t_mqrfr_total,t_cppRR_total,
                                 t_bmea_total,t_zippel_total,t_iratrecon_total,
                                 t_roots_total,
                                 t_genmono_total,t_construct_total,
                                 t_points_total,t_seqeval_total,
                                 t_normalize_total,t_NDSA_comp,
                                 stats_gc_count,stats_gc_Mbytes,
                                 solution_list,RUN_SYSTEM_ID,num_vars]]:
                else
                    printf("FAIL  n=%2d (mismatch)\n",n_test):
                    summary := [op(summary),
                                [n_test,"FAIL-mismatch",mrfi_calls,
                                 bb_calls_ndsa,bb_calls_mrfi,
                                 stats_terms_num,stats_terms_den,
                                 stats_deg_num,stats_deg_den,
                                 stats_common,stats_mMax,
                                 ffge_status,t_ffge,
                                 ffge_terms_num,ffge_terms_den,
                                 ffge_y_terms,ffge_det_terms,
                                 ffge_elim_swell,ffge_backsub_swell,
                                 t_ndsa_total,t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total,t_mqrfr_total,t_cppRR_total,
                                 t_bmea_total,t_zippel_total,t_iratrecon_total,
                                 t_roots_total,
                                 t_genmono_total,t_construct_total,
                                 t_points_total,t_seqeval_total,
                                 t_normalize_total,t_NDSA_comp,
                                 stats_gc_count,stats_gc_Mbytes,
                                 solution_list,RUN_SYSTEM_ID,num_vars]]:
                fi:
            catch:
                printf("PARTIAL  n=%2d  (ref solve threw)\n",n_test):
                for i from 1 to num_eqn do
                    printf("  x%d = %a\n", i, Final_rat_poly[i]):
                od:
                summary := [op(summary),
                            [n_test,"NO-VERIFY",mrfi_calls,
                             bb_calls_ndsa,bb_calls_mrfi,
                             stats_terms_num,stats_terms_den,
                             stats_deg_num,stats_deg_den,
                             stats_common,stats_mMax,
                                 ffge_status,t_ffge,
                                 ffge_terms_num,ffge_terms_den,
                                 ffge_y_terms,ffge_det_terms,
                                 ffge_elim_swell,ffge_backsub_swell,
                                 t_ndsa_total,t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total,t_mqrfr_total,t_cppRR_total,
                                 t_bmea_total,t_zippel_total,t_iratrecon_total,
                                 t_roots_total,
                                 t_genmono_total,t_construct_total,
                                 t_points_total,t_seqeval_total,
                                 t_normalize_total,t_NDSA_comp,
                                 stats_gc_count,stats_gc_Mbytes,
                                 solution_list,RUN_SYSTEM_ID,num_vars]]:
            end try:
        else
            printf("RECOVERED (unverified) n=%2d   BB-calls=%d\n",
                   n_test, mrfi_calls):
            for i from 1 to num_eqn do
                printf("  x%d = %a\n", i, Final_rat_poly[i]):
            od:
            summary := [op(summary),
                        [n_test,"UNVERIFIED",mrfi_calls,
                         bb_calls_ndsa,bb_calls_mrfi,
                         stats_terms_num,stats_terms_den,
                         stats_deg_num,stats_deg_den,
                         stats_common,stats_mMax,
                                 ffge_status,t_ffge,
                                 ffge_terms_num,ffge_terms_den,
                                 ffge_y_terms,ffge_det_terms,
                                 ffge_elim_swell,ffge_backsub_swell,
                                 t_ndsa_total,t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total,t_mqrfr_total,t_cppRR_total,
                                 t_bmea_total,t_zippel_total,t_iratrecon_total,
                                 t_roots_total,
                                 t_genmono_total,t_construct_total,
                                 t_points_total,t_seqeval_total,
                                 t_normalize_total,t_NDSA_comp,
                                 stats_gc_count,stats_gc_Mbytes,
                                 solution_list,RUN_SYSTEM_ID,num_vars]]:
        fi:
    else
        solution_list := []:
        printf("FAIL  n=%2d   %s\n",n_test,status):
        summary := [op(summary),
                    [n_test,status,mrfi_calls,
                     bb_calls_ndsa,bb_calls_mrfi,
                     stats_terms_num,stats_terms_den,
                     stats_deg_num,stats_deg_den,
                     stats_common,stats_mMax,
                                 ffge_status,t_ffge,
                                 ffge_terms_num,ffge_terms_den,
                                 ffge_y_terms,ffge_det_terms,
                                 ffge_elim_swell,ffge_backsub_swell,
                                 t_ndsa_total,t_mrfi_total,
                                 t_mrfi_wall,
                                 t_cppNewton_total,t_mqrfr_total,t_cppRR_total,
                                 t_bmea_total,t_zippel_total,t_iratrecon_total,
                                 t_roots_total,
                                 t_genmono_total,t_construct_total,
                                 t_points_total,t_seqeval_total,
                                 t_normalize_total,t_NDSA_comp,
                                 stats_gc_count,stats_gc_Mbytes,
                                 solution_list,RUN_SYSTEM_ID,num_vars]]:
    fi:
od:

(* 
print("=================================================================="):
print("  SUMMARY"):
print("=================================================================="):
printf("%4s  %-15s  %10s  %22s  %22s\n",
       "n", "status", "BB-calls",
       "MRFI num terms", "MRFI den terms"):
for entry in summary do
    printf("%4d  %-15s  %10d  %22a  %22a\n",
           entry[1], entry[2], entry[3],
           entry[6], entry[7]):
od:

if do_ffge then
    print(""):
    print("  Term count comparison: MRFI vs Lipson FFGE"):
    printf("%4s  %-22s  %-22s  %-22s  %-22s\n",
           "n",
           "MRFI num terms", "FFGE f terms",
           "MRFI den terms", "FFGE g terms"):
    for entry in summary do
        printf("%4d  %-22a  %-22a  %-22a  %-22a\n",
               entry[1], entry[6], entry[14], entry[7], entry[15]):
    od:
fi:
*)

report_path := "/cecm/home/mss59/Desktop/resMaple_MS/rat_recon_mpl/timings/FTR_sys_timing.txt":
# report_path := "/home/msokhi/Desktop/research_MS/rat_recon_mpl/timings/FTR_sys_Timing.txt":

# Label the report from the data themselves, not from the mutable selector.
REPORT_SYSTEM_ID := RUN_SYSTEM_ID:
if nops(summary) > 0 then
    REPORT_SYSTEM_ID := summary[1][39]:
    for entry in summary do
        if entry[39] <> REPORT_SYSTEM_ID then
            error "mixed system IDs in summary: expected %1 but found %2",
                  REPORT_SYSTEM_ID, entry[39]:
        fi:
    od:
fi:

fd := fopen(report_path, WRITE):
fprintf(fd, "============================================================\n"):
fprintf(fd, "  MRFI benchmark (fault-tolerant) -- system %s : %s\n",
        REPORT_SYSTEM_ID,
        ParametricSystems:-Name(REPORT_SYSTEM_ID)):
fprintf(fd, "  Univariate RR : HFTRFR (degrees + bad set) + DFTRFR\n"):
fprintf(fd, "  Prime p = %d\n", test_prime):
fprintf(fd, "  Range   n = %d .. %d\n", n_min, n_max):
fprintf(fd, "  Fault injection FAULT_ON = %a , error budget FAULT_E = %d\n",
        FAULT_ON, FAULT_E):
fprintf(fd, "============================================================\n\n"):

for entry in summary do
    t_wall  := entry[22]:
    t_comp  := entry[23]+entry[24]+entry[25]+entry[26]+entry[27]+entry[28]
              +entry[29]+entry[30]+entry[31]+entry[32]+entry[33]+entry[34]
              +entry[35]:
    t_over  := t_wall - t_comp:

    fprintf(fd, "  n = %d\n", entry[1]):
    fprintf(fd, "  System ID                        : %s\n", entry[39]):
    fprintf(fd, "  Parameter count                  : %d\n", entry[40]):
    fprintf(fd, "  MRFI Status                      : %s\n", entry[2]):
    fprintf(fd, "  Total BB calls                   : %d  (NDSA = %d, MRFI = %d)\n",
            entry[3], entry[4], entry[5]):
    fprintf(fd, "  BB calls in NDSA                 : %d\n", entry[4]):
    fprintf(fd, "  BB calls in MRFI                 : %d\n", entry[5]):
    fprintf(fd, "  Total NDSA time for BB calls (s) : %.9f\n", entry[20]):
    fprintf(fd, "  Total MRFI time for BB calls (s) : %.9f\n", entry[21]):
    fprintf(fd, "  Total BB calls time (NDSA+MRFI)  : %.9f\n", entry[20]+entry[21]):
    fprintf(fd, "  Total time (s)                   : %.9f\n", t_wall):
    fprintf(fd, "  Time for MRFI Individual component -\n"):
    fprintf(fd, "  Time -> cppNewtonInterp (s)  : %.9f\n", entry[23]):
    fprintf(fd, "  Time -> HFTRFR (s)           : %.9f\n", entry[24]):
    fprintf(fd, "  Time -> FTREval Interp+RR+EV : %.9f\n", entry[25]):
    fprintf(fd, "  Time -> BMEA (s)             : %.9f\n", entry[26]):
    fprintf(fd, "  Time -> Vandermonde (s)      : %.9f\n", entry[27]):
    fprintf(fd, "  Time -> Roots (s)            : %.9f\n", entry[29]):
    fprintf(fd, "  Time -> Gen. Monomials (s)   : %.9f\n", entry[30]):
    fprintf(fd, "  Time -> Final Poly. (s)      : %.9f\n", entry[31]):
    fprintf(fd, "  Time -> Get point AFL (s)    : %.9f\n", entry[32]):
    fprintf(fd, "  Time -> Eval+List (s)        : %.9f\n", entry[33]):
    fprintf(fd, "  Time -> Normalize Proc. (s)  : %.9f\n", entry[34]):
    fprintf(fd, "  Time -> MapRR (s)            : %.9f\n", entry[28]):
    fprintf(fd, "  Time -> Overhead & Other (s) : %.9f\n", t_over):
    fprintf(fd, "  Time -> NDSA Proc. (s)       : %.9f\n", entry[35]):
    fprintf(fd, "  GC -> collections            : %d\n", entry[36]):
    fprintf(fd, "  GC -> Mbytes allocated       : %.3f\n", entry[37]):
    fprintf(fd, "  Deg_num (per equation)   : %a\n", entry[8]):
    fprintf(fd, "  Deg_den (per equation)   : %a\n", entry[9]):
    fprintf(fd, "  Terms_num (per equation) : %a\n", entry[6]):
    fprintf(fd, "  Terms_den (per equation) : %a\n", entry[7]):
    fprintf(fd, "\n"):

    #  Fault tolerance overhead (per equation, 2E extra black-box points)
    
    fprintf(fd, "  Fault tolerance\n"):
    fprintf(fd, "  Error budget E               : %d\n", FAULT_E):
    fprintf(fd, "  Fault injection active       : %a\n", FAULT_ON):
    fprintf(fd, "  Extra BB pts / eqn (2E)      : %d\n", 2*FAULT_E):
    fprintf(fd, "\n"):
    if nops(entry[38]) > 0 then
        fprintf(fd, "  Recovered solution\n"):
        for i from 1 to nops(entry[38]) do
            fprintf(fd, "    x%d = %a\n", i, entry[38][i]):
        od:
    fi:
    fprintf(fd, "\n"):
    fprintf(fd, "  FFGE status              : %s\n", entry[12]):
    if entry[12] <> "SKIPPED" then
        fprintf(fd, "  FFGE total time (s)      : %.9f\n", entry[13]):
        fprintf(fd, "  FFGE f[i] terms          : %a\n", entry[14]):
        fprintf(fd, "  FFGE g[i] terms          : %a\n", entry[15]):
        fprintf(fd, "  FFGE y[i] terms (Pre GCD): %a\n", entry[16]):
        fprintf(fd, "  FFGE det(A) terms        : %a\n", entry[17]):
        fprintf(fd, "  FFGE Max Elim. step swell: %a  (numterms(num) before exact div)\n",
                entry[18]):
        fprintf(fd, "  FFGE max Back sub swell  : %a  (numterms(N[i]) before exact div)\n",
                entry[19]):
        fprintf(fd, "  Term count comparison (MRFI vs FFGE)\n"):
        fprintf(fd, "  Numerator terms MRFI     : %a\n", entry[6]):
        fprintf(fd, "  Numerator terms FFGE     : %a\n", entry[14]):
        fprintf(fd, "  Denominator terms MRFI   : %a\n", entry[7]):
        fprintf(fd, "  Denominator terms FFGE   : %a\n", entry[15]):
    fi:
    fprintf(fd, "\n"):
od:
fclose(fd):
printf("Wrote %s\n", report_path):
