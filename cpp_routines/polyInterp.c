/*
int NewtonInterp(LONG *X, LONG *Y, int n, LONG p) {    
// Interpolate f(X[i])=Y[i] for i=0,1,...,n-1.  
// The interpolating polynomial f is output Y.
     int d,i,j;
     // LONG prod,t,s,*V;
     // recint P = recip1(p);
     if( n<1 ) return -1;
     // V = Y;
     // Let f(z) = V_0 + V_1 (z-X_0) + ... + V_{n-1} (z-X_0) x ... x (z_X_{n-2})
     V[0] = Y[0];
     for (j=1; j<n; j++)
     {   s = V[0];
         prod = sub64s(X[j],X[0],p);
         for (i=1; i<j; i++)
         {   s = add64s(s,mulrec64(prod, V[i], P),p); // s += V[i] prod
             prod = mulrec64(prod,sub64s(X[j],X[i],p),P); // prod *= (X[j]-X[i)
         }
         if( prod==0 ) { printf("x coordinates must be distinct\n"); exit(1); }
         V[j] = mulrec64(sub64s(Y[j],s,p),modinv64s(prod,p),P);
     }
     for( d=n-1; d>=0 && V[d]==0; d-- ); // d = deg(f)
     for (i=1; i<=d; i++) // convert to standard basis using inplace Horner 
         for (j=d-i; j<=d-1; j++)
             V[j] = sub64s(V[j],mulrec64(X[d-i],V[j+1],P),p);
     return d;
}
*/ 