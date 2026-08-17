ParametricSystems := module()
    option package:

    export Generate,RawParams,ParamCount,DBSystem,Name,List,Catalog,DefaultSize:
    local catalog,order,defaults,path_stiffness,grid_laplacian,region_index,kron:

    region_index := proc(i::posint,r::posint)
        return ((i-1) mod r)+1:
    end proc:

    path_stiffness := proc(n::posint,diagshift::integer)
        local A,i:
        A := Matrix(n,n,fill=0):
        for i from 1 to n do
            A[i,i] := 2+diagshift:
            if i < n then
                A[i,i+1] := -1:
                A[i+1,i] := -1:
            fi:
        od:
        return A:
    end proc:

    grid_laplacian := proc(nv::posint,edges::list,wname)
        local L,e,u,v,ww,ground:
        L := Matrix(nv,nv,fill=0):
        for e from 1 to nops(edges) do
            u := edges[e][1]:
            v := edges[e][2]:
            ww := wname[e]:
            L[u,u] := L[u,u]+ww:
            L[v,v] := L[v,v]+ww:
            L[u,v] := L[u,v]-ww:
            L[v,u] := L[v,u]-ww:
        od:
        ground := nv:
        if nv=1 then
            return Matrix(0,0):
        fi:
        return L[1..ground-1,1..ground-1]:
    end proc:

    # Kept local so the database does not depend on a global with(LinearAlgebra).
    kron := proc(A::Matrix,B::Matrix)
        return LinearAlgebra:-KroneckerProduct(A,B):
    end proc:

    # ---- S: structured matrix families ------------------------------------

    catalog["S1"] := ["Symmetric Toeplitz",
        proc(n::posint)
            [Matrix(n,n,(i,j)->t[abs(i-j)]), Vector(n,fill=1)]
        end proc,
        proc(n::posint) [seq(t[i],i=0..n-1)] end proc]:

    catalog["S2"] := ["General Toeplitz",
        proc(n::posint)
            [Matrix(n,n,(i,j)->t[j-i]), Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(t[i],i=-(n-1)..n-1)] end proc]:

    catalog["S3"] := ["Hankel matrix",
        proc(n::posint)
            [Matrix(n,n,(i,j)->h[i+j-2]), Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(h[i],i=0..2*n-2)] end proc]:

    catalog["S4"] := ["Banded Toeplitz",
        proc(n::posint)
            local ww:
            ww := min(2,n-1):
            [Matrix(n,n,(i,j)->`if`(abs(j-i)<=ww,t[j-i],0)), Vector(n,i->i)]
        end proc,
        proc(n::posint)
            local ww:
            ww := min(2,n-1):
            [seq(t[i],i=-ww..ww)]
        end proc]:

    catalog["S5"] := ["Block Toeplitz with Toeplitz blocks (BTTB)",
        proc(q::posint)
            local N:
            N := q^2:
            [Matrix(N,N,
                    (i,j)->tau[floor((j-1)/q)-floor((i-1)/q),
                                ((j-1) mod q)-((i-1) mod q)]),
             Vector(N,i->i)]
        end proc,
        proc(q::posint)
            [seq(seq(tau[a,b],b=-(q-1)..q-1),a=-(q-1)..q-1)]
        end proc]:

    catalog["S6"] := ["Generic tridiagonal / continuant family",
        proc(n::posint)
            [Matrix(n,n,(i,j)->piecewise(i=j,d[i],j=i+1,u[i],i=j+1,l[j],0)),
             Vector(n,i->i)]
        end proc,
        proc(n::posint)
            [seq(d[i],i=1..n),seq(u[i],i=1..n-1),seq(l[i],i=1..n-1)]
        end proc]:

    # ---- R: algebraic / resultant families --------------------------------

    catalog["R1"] := ["Sylvester resultant matrix",
        proc(n::posint)
            local A,row,k,col:
            # Balanced case r=s=n, so N=2n.
            A := Matrix(2*n,2*n,fill=0):
            for row from 1 to n do
                for k from 0 to n do
                    col := row+k:
                    A[row,col] := a[n-k]:
                    A[n+row,col] := bcoef[n-k]:
                od:
            od:
            [A,Vector(2*n,i->i)]
        end proc,
        proc(n::posint)
            [seq(a[i],i=0..n),seq(bcoef[i],i=0..n)]
        end proc]:

    catalog["R2"] := ["Bezout matrix of two generic polynomials",
        proc(q::posint)
            local Bm,p,qq,k,ii,jj,cc:
            Bm := Matrix(q,q,fill=0):
            # For p>qq,
            # (x^p y^qq-y^p x^qq)/(x-y)
            # = sum(k=0..p-qq-1) x^(p-1-k)y^(qq+k).
            for p from 1 to q do
                for qq from 0 to p-1 do
                    cc := a[p]*bcoef[qq]-a[qq]*bcoef[p]:
                    for k from 0 to p-qq-1 do
                        ii := p-k:
                        jj := qq+k+1:
                        Bm[ii,jj] := Bm[ii,jj]+cc:
                    od:
                od:
            od:
            [Bm,Vector(q,i->i)]
        end proc,
        proc(q::posint)
            [seq(a[i],i=0..q),seq(bcoef[i],i=0..q)]
        end proc]:

    catalog["R3"] := ["Companion resolvent",
        proc(n::posint)
            local C,A,i:
            C := Matrix(n,n,fill=0):
            for i from 1 to n-1 do
                C[i+1,i] := 1:
            od:
            for i from 1 to n do
                C[i,n] := -a[i-1]:
            od:
            A := s*LinearAlgebra:-IdentityMatrix(n)-C:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [s,seq(a[i],i=0..n-1)] end proc]:

    # ---- P: application-motivated scalable families ------------------------

    catalog["P1"] := ["Grounded weighted q x q grid Laplacian",
        proc(q::posint)
            local N,L,node,i,j,ww:
            N := q^2:
            L := Matrix(N,N,fill=0):
            node := (i,j)->(i-1)*q+j:
            for i from 1 to q do
                for j from 1 to q-1 do
                    ww := wh[i,j]:
                    L[node(i,j),node(i,j)] := L[node(i,j),node(i,j)]+ww:
                    L[node(i,j+1),node(i,j+1)] := L[node(i,j+1),node(i,j+1)]+ww:
                    L[node(i,j),node(i,j+1)] := L[node(i,j),node(i,j+1)]-ww:
                    L[node(i,j+1),node(i,j)] := L[node(i,j+1),node(i,j)]-ww:
                od:
            od:
            for i from 1 to q-1 do
                for j from 1 to q do
                    ww := wv[i,j]:
                    L[node(i,j),node(i,j)] := L[node(i,j),node(i,j)]+ww:
                    L[node(i+1,j),node(i+1,j)] := L[node(i+1,j),node(i+1,j)]+ww:
                    L[node(i,j),node(i+1,j)] := L[node(i,j),node(i+1,j)]-ww:
                    L[node(i+1,j),node(i,j)] := L[node(i+1,j),node(i,j)]-ww:
                od:
            od:
            [L[1..N-1,1..N-1],Vector(N-1,i->piecewise(i=1,1,0))]
        end proc,
        proc(q::posint)
            [seq(seq(wh[i,j],j=1..q-1),i=1..q),
             seq(seq(wv[i,j],j=1..q),i=1..q-1)]
        end proc]:

    catalog["P2"] := ["2-D convection-diffusion-reaction finite difference",
        proc(q::posint)
            local N,A,node,i,j,k,src_i,src_j:
            N := q^2:
            A := Matrix(N,N,fill=0):
            node := (i,j)->(i-1)*q+j:
            for i from 1 to q do
                for j from 1 to q do
                    k := node(i,j):
                    A[k,k] := 8+2*gamma:
                    if j<q then A[k,node(i,j+1)] := -2+alpha fi:
                    if j>1 then A[k,node(i,j-1)] := -2-alpha fi:
                    if i<q then A[k,node(i+1,j)] := -2+beta fi:
                    if i>1 then A[k,node(i-1,j)] := -2-beta fi:
                od:
            od:
            src_i := iquo(q+1,2): if src_i<1 then src_i:=1 fi:
            src_j := src_i:
            [A,Vector(N,k->piecewise(k=node(src_i,src_j),1,0))]
        end proc,
        proc(q::posint) [alpha,beta,gamma] end proc]:

    catalog["P3"] := ["RC ladder nodal admittance matrix",
        proc(n::posint)
            local A,i,y1,y2:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                y1 := g[i]+s*c[i]:
                y2 := g[i+1]+s*c[i+1]:
                A[i,i] := y1+y2:
                if i<n then
                    A[i,i+1] := -y2:
                    A[i+1,i] := -y2:
                fi:
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) [s,seq(g[i],i=1..n+1),seq(c[i],i=1..n+1)] end proc]:

    catalog["P4"] := ["Mass-spring chain dynamic stiffness",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := k[i]+k[i+1]-dinert[i]:
                if i<n then
                    A[i,i+1] := -k[i+1]:
                    A[i+1,i] := -k[i+1]:
                fi:
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) [seq(k[i],i=1..n+1),seq(dinert[i],i=1..n)] end proc]:

    catalog["P5"] := ["Birth-death Markov-chain reachability",
        proc(n::posint)
            local A,rhs,i:
            A := Matrix(n,n,fill=0):
            rhs := Vector(n,fill=0):
            for i from 1 to n do
                A[i,i] := p[i]+qrate[i]:
                if i<n then A[i,i+1] := -p[i] else rhs[i] := p[i] fi:
                if i>1 then A[i,i-1] := -qrate[i] fi:
            od:
            [A,rhs]
        end proc,
        proc(n::posint) [seq(p[i],i=1..n),seq(qrate[i],i=1..n)] end proc]:

    catalog["P6"] := ["General grounded resistor / conductance network",
        proc(n::posint)
            local edges,i,A:
            edges := [seq([i,i+1],i=1..n),seq([i,i+2],i=1..n-1)]:
            A := grid_laplacian(n+1,edges,gnet):
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) [seq(gnet[i],i=1..2*n-1)] end proc]:

    catalog["P7"] := ["DC power-flow network",
        proc(n::posint)
            local edges,i,A:
            edges := [seq([i,i+1],i=1..n),seq([i,i+2],i=1..n-1)]:
            A := grid_laplacian(n+1,edges,bsus):
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(bsus[i],i=1..2*n-1)] end proc]:

    catalog["P8"] := ["1-D variable-conductivity diffusion / heat equation",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := kap[i]+kap[i+1]:
                if i>1 then A[i,i-1] := -kap[i] fi:
                if i<n then A[i,i+1] := -kap[i+1] fi:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(kap[i],i=1..n+1)] end proc]:

    catalog["P9"] := ["2-D anisotropic variable-coefficient diffusion",
        proc(q::posint)
            local N,A,node,i,j,k,west,east,south,north:
            N := q^2:
            A := Matrix(N,N,fill=0):
            node := (i,j)->(i-1)*q+j:
            for i from 1 to q do
                for j from 1 to q do
                    k := node(i,j):
                    west  := ax[i,j-1]:
                    east  := ax[i,j]:
                    south := ay[i-1,j]:
                    north := ay[i,j]:
                    A[k,k] := west+east+south+north:
                    if j>1 then A[k,node(i,j-1)] := -west fi:
                    if j<q then A[k,node(i,j+1)] := -east fi:
                    if i>1 then A[k,node(i-1,j)] := -south fi:
                    if i<q then A[k,node(i+1,j)] := -north fi:
                od:
            od:
            [A,Vector(N,i->i)]
        end proc,
        proc(q::posint)
            [seq(seq(ax[i,j],j=0..q),i=1..q),
             seq(seq(ay[i,j],j=1..q),i=0..q)]
        end proc]:

    catalog["P10"] := ["Thermal model with parametric Robin / film coefficients",
        proc(n::posint)
            local A,r,i,j:
            r := min(3,n):
            A := path_stiffness(n,1):
            for i from 1 to n do
                j := region_index(i,r):
                A[i,i] := A[i,i]+hfilm[j]:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[seq(hfilm[j],j=1..r)] end proc]:

    catalog["P11"] := ["Transient thermal descriptor / Laplace-domain heat model",
        proc(n::posint)
            local A,r,i,j:
            r := min(3,n):
            A := path_stiffness(n,1):
            for i from 1 to n do
                j := region_index(i,r):
                A[i,i] := A[i,i]+kmat[j]+s*(1+((i-1) mod 3)):
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[s,seq(kmat[j],j=1..r)] end proc]:

    catalog["P12"] := ["Darcy / groundwater / single-phase reservoir pressure",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := trans[i]+trans[i+1]+well[i]:
                if i>1 then A[i,i-1] := -trans[i] fi:
                if i<n then A[i,i+1] := -trans[i+1] fi:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(trans[i],i=1..n+1),seq(well[i],i=1..n)] end proc]:

    catalog["P13"] := ["Variable-viscosity Stokes mixed finite-element system",
        proc(n::posint)
            local np,N,A,r,i,j,row:
            if n<2 then error "P13 requires n >= 2" fi:
            np := n-1: N := n+np: r := min(4,n):
            A := Matrix(N,N,fill=0):
            for i from 1 to n do
                A[i,i] := nu[region_index(i,r)]:
            od:
            for row from 1 to np do
                A[n+row,row] := -1: A[n+row,row+1] := 1:
                A[row,n+row] := -1: A[row+1,n+row] := 1:
            od:
            [A,Vector(N,i->piecewise(i<=n,i,0))]
        end proc,
        proc(n::posint) local r:if n<2 then error "P13 requires n >= 2" fi:r:=min(4,n):[seq(nu[j],j=1..r)] end proc]:

    catalog["P14"] := ["Mixed Darcy saddle-point system",
        proc(n::posint)
            local np,N,A,r,i,j,row:
            if n<2 then error "P14 requires n >= 2" fi:
            np := n-1: N := n+np: r := min(4,n):
            A := Matrix(N,N,fill=0):
            for i from 1 to n do
                A[i,i] := rho[region_index(i,r)]:
            od:
            for row from 1 to np do
                A[n+row,row] := -1: A[n+row,row+1] := 1:
                A[row,n+row] := -1: A[row+1,n+row] := 1:
            od:
            [A,Vector(N,i->piecewise(i<=n,i,1))]
        end proc,
        proc(n::posint) local r:if n<2 then error "P14 requires n >= 2" fi:r:=min(4,n):[seq(rho[j],j=1..r)] end proc]:

    catalog["P15"] := ["Helmholtz / acoustic frequency-response system",
        proc(n::posint)
            local A,r,i,j:
            r := min(3,n): A := path_stiffness(n,0):
            for i from 1 to n do
                j := region_index(i,r):
                A[i,i] := A[i,i]+khelm[j]-zfreq*(1+((i-1) mod 2)):
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[zfreq,seq(khelm[j],j=1..r)] end proc]:

    catalog["P16"] := ["Static truss / linear-elasticity stiffness system",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := thetaE[i]+thetaE[i+1]:
                if i>1 then A[i,i-1] := -thetaE[i] fi:
                if i<n then A[i,i+1] := -thetaE[i+1] fi:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(thetaE[i],i=1..n+1)] end proc]:

    catalog["P17"] := ["Structural FEM frequency response",
        proc(n::posint)
            local A,r,i:
            r := min(3,n): A := path_stiffness(n,1):
            for i from 1 to n do
                A[i,i] := A[i,i]+thetaS[region_index(i,r)]-zfreq*(1+((i-1) mod 3)):
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[zfreq,seq(thetaS[j],j=1..r)] end proc]:

    catalog["P18"] := ["Damped structural dynamics affine dynamic stiffness",
        proc(n::posint)
            local A,r,i:
            r := min(3,n): A := path_stiffness(n,1):
            for i from 1 to n do
                A[i,i] := A[i,i]+thetaS[region_index(i,r)]+2*dshift+zshift*(1+((i-1) mod 3)):
                if i<n then
                    A[i,i+1] := A[i,i+1]-dshift:
                    A[i+1,i] := A[i+1,i]-dshift:
                fi:
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[seq(thetaS[j],j=1..r),dshift,zshift] end proc]:

    catalog["P19"] := ["RLC modified-nodal-analysis network",
        proc(n::posint)
            local A,N,i:
            N := 2*n: A := Matrix(N,N,fill=0):
            for i from 1 to n do
                A[i,i] := gRLC[i]+etaRLC[i]:
                A[n+i,n+i] := -ellRLC[i]:
                A[i,n+i] := 1: A[n+i,i] := 1:
                if i>1 then
                    A[i,n+i-1] := -1:
                    A[n+i-1,i] := -1:
                fi:
            od:
            [A,Vector(N,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint)
            [seq(gRLC[i],i=1..n),seq(etaRLC[i],i=1..n),seq(ellRLC[i],i=1..n)]
        end proc]:

    catalog["P20"] := ["Eddy-current / low-frequency electromagnetic FEM",
        proc(n::posint)
            local A,r1,r2,i:
            r1 := min(3,n): r2 := min(3,n): A := path_stiffness(n,1):
            for i from 1 to n do
                A[i,i] := A[i,i]+reluct[region_index(i,r1)]+etaEM[region_index(i,r2)]:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) local r1,r2:r1:=min(3,n):r2:=min(3,n):[seq(reluct[j],j=1..r1),seq(etaEM[j],j=1..r2)] end proc]:

    catalog["P21"] := ["Fixed-source discrete-ordinates / radiative transport",
        proc(n::posint)
            local A,N,c,ip,im:
            N := 2*n: A := Matrix(N,N,fill=0):
            for c from 1 to n do
                ip := c: im := n+c:
                A[ip,ip] := 1+sigT[c]:
                A[im,im] := 1+sigT[c]:
                A[ip,im] := -sigS[c]:
                A[im,ip] := -sigS[c]:
                if c>1 then A[ip,ip-1] := -1 fi:
                if c<n then A[im,im+1] := -1 fi:
            od:
            [A,Vector(N,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) [seq(sigT[i],i=1..n),seq(sigS[i],i=1..n)] end proc]:

    catalog["P22"] := ["Multigroup neutron-diffusion fixed-source system",
        proc(n::posint)
            local K,A,N,i,j:
            K := path_stiffness(n,0): N := 2*n: A := Matrix(N,N,fill=0):
            for i from 1 to n do
                for j from 1 to n do
                    A[i,j] := Dg[1]*K[i,j]:
                    A[n+i,n+j] := Dg[2]*K[i,j]:
                od:
                A[i,i] := A[i,i]+remg[1]:
                A[n+i,n+i] := A[n+i,n+i]+remg[2]:
                A[i,n+i] := -scat21:
                A[n+i,i] := -scat12:
            od:
            [A,Vector(N,i->i)]
        end proc,
        proc(n::posint) [Dg[1],Dg[2],remg[1],remg[2],scat12,scat21] end proc]:

    catalog["P23"] := ["General absorbing Markov-chain reachability",
        proc(n::posint)
            local A,rhs,i:
            A := LinearAlgebra:-IdentityMatrix(n): rhs := Vector(n,fill=0):
            for i from 1 to n do
                if i<n then A[i,i+1] := -pr[i] fi:
                if i>1 then A[i,i-1] := -pl[i] fi:
                rhs[i] := pa[i]:
            od:
            [A,rhs]
        end proc,
        proc(n::posint) [seq(pr[i],i=1..n-1),seq(pl[i],i=2..n),seq(pa[i],i=1..n)] end proc]:

    catalog["P24"] := ["Continuous-time Markov mean first-passage time",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := absorb[i]:
                if i<n then A[i,i] := A[i,i]+lam[i]: A[i,i+1] := -lam[i] fi:
                if i>1 then A[i,i] := A[i,i]+mu[i]: A[i,i-1] := -mu[i] fi:
            od:
            [A,Vector(n,fill=1)]
        end proc,
        proc(n::posint) [seq(lam[i],i=1..n-1),seq(mu[i],i=2..n),seq(absorb[i],i=1..n)] end proc]:

    catalog["P25"] := ["CTMC stationary distribution with normalization",
        proc(n::posint)
            local Q,A,rhs,i,jf,jb:
            if n<2 then error "P25 requires n >= 2" fi:
            Q := Matrix(n,n,fill=0):
            for i from 1 to n do
                jf := (i mod n)+1:
                jb := ((i-2) mod n)+1:
                Q[i,jf] := Q[i,jf]+kf[i]:
                Q[i,jb] := Q[i,jb]+kb[i]:
                Q[i,i] := Q[i,i]-kf[i]-kb[i]:
            od:
            A := LinearAlgebra:-Transpose(Q):
            for i from 1 to n do A[n,i] := 1 od:
            rhs := Vector(n,fill=0): rhs[n] := 1:
            [A,rhs]
        end proc,
        proc(n::posint) if n<2 then error "P25 requires n >= 2" fi:[seq(kf[i],i=1..n),seq(kb[i],i=1..n)] end proc]:

    catalog["P26"] := ["PageRank / personalized PageRank linear system",
        proc(n::posint)
            local Pm,A,rhs,i,j:
            Pm := Matrix(n,n,fill=0):
            for i from 1 to n do
                j := (i mod n)+1:
                Pm[i,j] := 1:
            od:
            A := LinearAlgebra:-IdentityMatrix(n)-pagerank_alpha*LinearAlgebra:-Transpose(Pm):
            rhs := Vector(n,i->piecewise(i=1,1-pagerank_alpha,0)):
            [A,rhs]
        end proc,
        proc(n::posint) [pagerank_alpha] end proc]:

    catalog["P27"] := ["Discounted Markov reward / fixed-policy value system",
        proc(n::posint)
            local Pm,A,i,j:
            Pm := Matrix(n,n,fill=0):
            for i from 1 to n do
                j := (i mod n)+1:
                Pm[i,j] := 1:
            od:
            A := LinearAlgebra:-IdentityMatrix(n)-discount_gamma*Pm:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [discount_gamma] end proc]:

    catalog["P28"] := ["First-order reaction / compartment / pharmacokinetic network",
        proc(n::posint)
            local A,i,j:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                j := (i mod n)+1:
                A[i,i] := sink[i]+ktr[i]:
                A[j,i] := A[j,i]-ktr[i]:
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) [seq(ktr[i],i=1..n),seq(sink[i],i=1..n)] end proc]:

    catalog["P29"] := ["Leontief input-output production network",
        proc(n::posint)
            local A,i:
            A := LinearAlgebra:-IdentityMatrix(n):
            for i from 1 to n do
                A[i,i] := A[i,i]-cdiag[i]:
                if i<n then
                    A[i,i+1] := -cup[i]:
                    A[i+1,i] := -clow[i]:
                fi:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(cdiag[i],i=1..n),seq(cup[i],i=1..n-1),seq(clow[i],i=1..n-1)] end proc]:

    catalog["P30"] := ["Weighted least squares / sensor fusion normal equations",
        proc(n::posint)
            local m,A,rhs,obs,i,j,idx1,idx2,zobs:
            m := 2*n-1:
            A := Matrix(n,n,fill=0): rhs := Vector(n,fill=0):
            # First n observations have h=e_i.
            for obs from 1 to n do
                A[obs,obs] := A[obs,obs]+wobs[obs]:
                rhs[obs] := rhs[obs]+wobs[obs]*obs:
            od:
            # Next n-1 observations have h=e_i+e_{i+1}.
            for i from 1 to n-1 do
                obs := n+i: zobs := n+i:
                A[i,i] := A[i,i]+wobs[obs]:
                A[i+1,i+1] := A[i+1,i+1]+wobs[obs]:
                A[i,i+1] := A[i,i+1]+wobs[obs]:
                A[i+1,i] := A[i+1,i]+wobs[obs]:
                rhs[i] := rhs[i]+wobs[obs]*zobs:
                rhs[i+1] := rhs[i+1]+wobs[obs]*zobs:
            od:
            for i from 1 to n do A[i,i] := A[i,i]+reglambda od:
            [A,rhs]
        end proc,
        proc(n::posint) [seq(wobs[i],i=1..2*n-1),reglambda] end proc]:

    catalog["P31"] := ["Equality-constrained quadratic-program KKT system",
        proc(n::posint)
            local m,N,A,row,i:
            if n<2 then error "P31 requires n >= 2" fi:
            m := max(1,iquo(n,2)): N := n+m: A := Matrix(N,N,fill=0):
            for i from 1 to n do A[i,i] := 2+thetaQP[i] od:
            for row from 1 to m do
                A[n+row,row] := 1: A[n+row,row+1] := -1:
                A[row,n+row] := 1: A[row+1,n+row] := -1:
            od:
            [A,Vector(N,i->piecewise(i<=n,-i,1))]
        end proc,
        proc(n::posint) if n<2 then error "P31 requires n >= 2" fi:[seq(thetaQP[i],i=1..n)] end proc]:

    catalog["P32"] := ["Interior-point / regularized normal-equation family",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := dbar[i]+dbar[i+1]+rhoReg:
                if i<n then
                    A[i,i+1] := -dbar[i+1]:
                    A[i+1,i] := -dbar[i+1]:
                fi:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [seq(dbar[i],i=1..n+1),rhoReg] end proc]:

    catalog["P33"] := ["Continuous-time Lyapunov equation (vectorized)",
        proc(n::posint)
            local Am,Id,K,r,i,j,Qv:
            r := min(3,n): Am := path_stiffness(n,1):
            for i from 1 to n do Am[i,i] := Am[i,i]+thetaL[region_index(i,r)] od:
            Id := LinearAlgebra:-IdentityMatrix(n):
            K := kron(Id,LinearAlgebra:-Transpose(Am))
                 +kron(LinearAlgebra:-Transpose(Am),Id):
            Qv := Vector(n^2,idx->piecewise(iquo(idx-1,n)+1=((idx-1) mod n)+1,-1,0)):
            [K,Qv]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[seq(thetaL[j],j=1..r)] end proc]:

    catalog["P34"] := ["Sylvester matrix equation (vectorized)",
        proc(n::posint)
            local AL,BR,Id,K,r1,r2,i:
            r1 := min(2,n): r2 := min(2,n):
            AL := path_stiffness(n,1): BR := path_stiffness(n,2):
            for i from 1 to n do
                AL[i,i] := AL[i,i]+thetaA[region_index(i,r1)]:
                BR[i,i] := BR[i,i]+muB[region_index(i,r2)]:
            od:
            Id := LinearAlgebra:-IdentityMatrix(n):
            K := kron(Id,AL)+kron(LinearAlgebra:-Transpose(BR),Id):
            [K,Vector(n^2,fill=1)]
        end proc,
        proc(n::posint) local r1,r2:r1:=min(2,n):r2:=min(2,n):[seq(thetaA[j],j=1..r1),seq(muB[j],j=1..r2)] end proc]:

    catalog["P35"] := ["Descriptor state-space transfer resolvent / parametric MOR",
        proc(n::posint)
            local A,r,i:
            r := min(3,n): A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := sMOR*(1+((i-1) mod 3))+2-thetaMOR[region_index(i,r)]:
                if i<n then
                    A[i,i+1] := -1:
                    A[i+1,i] := -1:
                fi:
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[sMOR,seq(thetaMOR[j],j=1..r)] end proc]:

    catalog["P36"] := ["Gaussian covariance-component / kriging solve",
        proc(n::posint)
            local A,r,i,j,k:
            r := min(3,n): A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := sigma2:
                for j from 1 to n do
                    for k from 1 to r do
                        A[i,j] := A[i,j]+thetaCov[k]*(i+k)*(j+k):
                    od:
                od:
            od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) local r:r:=min(3,n):[sigma2,seq(thetaCov[k],k=1..r)] end proc]:

    catalog["P37"] := ["Spline smoothing / penalized regression system",
        proc(n::posint)
            local A,rhs,i,obs,zobs:
            A := Matrix(n,n,fill=0): rhs := Vector(n,fill=0):
            # B is fixed upper bidiagonal: row i = e_i+e_{i+1} (last row e_n).
            for obs from 1 to n do
                zobs := obs:
                A[obs,obs] := A[obs,obs]+wspline[obs]:
                rhs[obs] := rhs[obs]+wspline[obs]*zobs:
                if obs<n then
                    A[obs+1,obs+1] := A[obs+1,obs+1]+wspline[obs]:
                    A[obs,obs+1] := A[obs,obs+1]+wspline[obs]:
                    A[obs+1,obs] := A[obs+1,obs]+wspline[obs]:
                    rhs[obs+1] := rhs[obs+1]+wspline[obs]*zobs:
                fi:
            od:
            # P = 1-D second-difference penalty.
            for i from 1 to n do
                A[i,i] := A[i,i]+2*lambdaSpline:
                if i<n then
                    A[i,i+1] := A[i,i+1]-lambdaSpline:
                    A[i+1,i] := A[i+1,i]-lambdaSpline:
                fi:
            od:
            [A,rhs]
        end proc,
        proc(n::posint) [seq(wspline[i],i=1..n),lambdaSpline] end proc]:

    catalog["P38"] := ["Tight-binding / discretized Schrodinger resolvent",
        proc(n::posint)
            local A,i:
            A := Matrix(n,n,fill=0):
            for i from 1 to n do
                A[i,i] := 2+vpot[i]-energyE:
                if i<n then
                    A[i,i+1] := -1:
                    A[i+1,i] := -1:
                fi:
            od:
            [A,Vector(n,i->piecewise(i=1,1,0))]
        end proc,
        proc(n::posint) [energyE,seq(vpot[i],i=1..n)] end proc]:

    catalog["P39"] := ["Age-structured / Leslie population resolvent",
        proc(n::posint)
            local A,i:
            A := zLeslie*LinearAlgebra:-IdentityMatrix(n):
            for i from 1 to n do A[1,i] := A[1,i]-fert[i] od:
            for i from 1 to n-1 do A[i+1,i] := A[i+1,i]-surv[i] od:
            [A,Vector(n,i->i)]
        end proc,
        proc(n::posint) [zLeslie,seq(fert[i],i=1..n),seq(surv[i],i=1..n-1)] end proc]:

    catalog["P40"] := ["Quasi-birth-death queueing reachability",
        proc(L::posint)
            local r,N,A,lev,i,j,row,col:
            r := 2: N := L*r: A := Matrix(N,N,fill=0):
            for lev from 1 to L do
                for i from 1 to r do
                    for j from 1 to r do
                        row := (lev-1)*r+i:
                        col := (lev-1)*r+j:
                        A[row,col] := piecewise(i=j,1,0)-a0q[i,j]:
                        if lev>1 then
                            col := (lev-2)*r+j:
                            A[row,col] := -aminus[i,j]:
                        fi:
                        if lev<L then
                            col := lev*r+j:
                            A[row,col] := -aplus[i,j]:
                        fi:
                    od:
                od:
            od:
            [A,Vector(N,i->piecewise(i=1,1,0))]
        end proc,
        proc(L::posint)
            [seq(seq(aminus[i,j],j=1..2),i=1..2),
             seq(seq(a0q[i,j],j=1..2),i=1..2),
             seq(seq(aplus[i,j],j=1..2),i=1..2)]
        end proc]:

    # ---- order and suggested first sizes ----------------------------------

    order := ["S1","S2","S3","S4","S5","S6",
              "R1","R2","R3",
              "P1","P2","P3","P4","P5","P6","P7","P8","P9","P10",
              "P11","P12","P13","P14","P15","P16","P17","P18","P19","P20",
              "P21","P22","P23","P24","P25","P26","P27","P28","P29","P30",
              "P31","P32","P33","P34","P35","P36","P37","P38","P39","P40"]:

    defaults := table([
        "S1"=8,"S2"=6,"S3"=6,"S4"=50,"S5"=4,"S6"=20,
        "R1"=8,"R2"=6,"R3"=50,
        "P1"=4,"P2"=8,"P3"=20,"P4"=20,"P5"=20
    ]):

    # ---- public API --------------------------------------------------------

    Generate := proc(id::string,n::posint)
        if not assigned(catalog[id]) then
            error "Unknown system ID %1. Run ParametricSystems:-List().",id:
        fi:
        return catalog[id][2](n):
    end proc:

    RawParams := proc(id::string,n::posint)
        if not assigned(catalog[id]) then
            error "Unknown system ID %1. Run ParametricSystems:-List().",id:
        fi:
        return catalog[id][3](n):
    end proc:

    ParamCount := proc(id::string,n::posint)
        return numelems(RawParams(id,n)):
    end proc:

    Name := proc(id::string)
        if not assigned(catalog[id]) then
            error "Unknown system ID %1. Run ParametricSystems:-List().",id:
        fi:
        return catalog[id][1]:
    end proc:

    DefaultSize := proc(id::string)
        if assigned(defaults[id]) then return defaults[id] fi:
        return 10:
    end proc:

    DBSystem := proc(id::string,n::posint)
        local res,A,b,nn,Sys,pars0,newp,remap,Vars,i,j:
        global DB_param_map,DB_xstar:

        res := Generate(id,n):
        A := res[1]: b := res[2]:
        nn := LinearAlgebra:-Dimension(b):
        Vars := [seq(x||i,i=1..nn)]:

        pars0 := RawParams(id,n):
        newp := [seq(y||i,i=1..numelems(pars0))]:
        remap := zip(`=`,pars0,newp):

        Sys := [seq(add(A[i,j]*Vars[j],j=1..nn)=b[i],i=1..nn)]:
        Sys := subs(remap,Sys):
        DB_param_map := remap:

        if numelems(res)>=3 then
            DB_xstar := subs(remap,convert(res[3],list)):
        else
            DB_xstar := 'DB_xstar':
        fi:

        return Sys,Vars,newp,numelems(newp),nn:
    end proc:

    List := proc()
        local id:
        for id in order do
            printf("%-4s  %s\n",id,catalog[id][1]):
        od:
        NULL:
    end proc:

    Catalog := proc()
        return eval(catalog):
    end proc:

end module:
