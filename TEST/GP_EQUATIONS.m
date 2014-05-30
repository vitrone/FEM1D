classdef GP_EQUATIONS

    methods (Static)
        function out = ABC1a_CQ1(Data)
            % Solving GP equation using FEM
            % iu_t+u_xx+xg(t)u+chi*|u|^2u=0
            % with artificial boundary conditions ABC_1.
            % ABC_1^(*) -- ABC_1 discretized using convolution quadrature (CQ).
            % Argument structure is 
            % Data
            %   |----fp  --->  p,     Highest degree of the polynomials
            %   |----fN  --->  N,     Number of partitions defining finite elements
            %   |----fx  --->  x,     Grid points
            %   |----fFM --->  FM,    Forward Transform matrix
            %   |----fIM --->  IM,    Inverse transform matrix
            %   |----fFMi--->  FMi,   Forward transform for anti-aliasing
            %   |----fIMi--->  IMi,   Inverse transform for anti-aliasing
            %   |----fP  --->  P,     Number of LGL points for Numerical integration
            %   |----fxri--->  xri,   LGL points for Numerical integration
            %   |----fquadWi---> quadWi,    Quadrature weights for LGL points
            %   |----ftol--->  tol,   tolerance
            %   |----fdt --->  dt,    Time step
            %   |----fNt --->  Nt,    Number of time steps
            %   |----fNth--->  Nth,   Threshold iteration
            %   |----fMp --->  Mp,    Order of the Pade approximants
            %   |----fx_r--->  x_r,   Left boundary
            %   |----fx_l--->  x_l,   Right boundary
            %   |----fAL --->  AL,    Value of the field on the left exterior domain
            %   |----fAR --->  AR,    Value of the field on the right exterior domain
            %   |----fKL --->  KL,    Wave vector for the left exterior domain
            %   |----fKR --->  KR,    Wave vector for the right exterior domain
            %   |----fchi--->  chi,   Nonlinear coefficient (refer to the evolution equation)
            %   |----fu0 --->  u0,    Initial data
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %
            % Degree of the Legendre polynomials involved --> p
            % tol = 1e-9;
            % [xr quadW IM FM] = LGLdataLT2( p, tol);
            % xr: LGL points
            % FM: matrix for forward transform
            % IM: matrix for inverse transform

            % Anti aliasing for the nonlinear term
            % P = 2*p;% For the purpose of numerical integration
            % Degree of the Legendre polynomials involved --> P-1
            % [xri quadWi] = LGLdataLT1( P, tol);
            % FMi = LGLdataFM( p, P, xri); 
            % IMi = LGLdataIM( P, p, xri);
            %
            % Pade method for fractional derivativesm requires 
            % the mex function PF_Pade.
            % 
            % Date : Apr 2013
            %
            % Updates
            %
            % Apr 2013: added this function as a class method
            %%
            FM      = Data.fFM;
            IM      = Data.fIM;
            FMi     = Data.fFMi;
            IMi     = Data.fIMi;
            quadWi  = Data.fquadWi;
            xri     = Data.fxri;
            tol     = Data.ftol;
            N       = Data.fN;
            p       = Data.fp;
            P       = Data.fP;
            dt      = Data.fdt;
            Nt      = Data.fNt;
            Nth     = Data.fNth;
            x_r     = Data.fx_r;
            x_l     = Data.fx_l;
            AL      = Data.fAL;
            AR      = Data.fAR;
            KL      = Data.fKL;
            KR      = Data.fKR;
            chi     = Data.fchi;
            x       = Data.fx;
            V       = str2func(Data.hpot);
            XiL     = str2func(Data.hXiL);
            XiR     = str2func(Data.hXiR);
            dXi     = str2func(Data.hdXi);
            profile = str2func(Data.hfunc);

            I  =  sqrt(-1);
            dx =  (x_r-x_l)/N;
            J  =  dx/2;
            r  =  2/dt;
             
            ph  =   exp(-I*pi/4);
            iph =   exp(I*pi/4);
            %% Convolution Quadrature
            G1o2  =  r^(1/2)*gfilters(1/2,Nt);
            r1o2  =  r^(1/2);
            Gm1o2 =  r^(-1/2)*gfilters(-1/2,Nt);
            rm1o2 =  r^(-1/2);
            %%
            Q = mxfem1d_quadM( quadWi, IMi);
            x1 = ref2mesh (xri, N, [x_l x_r]);
            % mass matrix
            pot = ones(size(x1));
            
            LGMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGMM = transpose(LGMM) - diag(diag(LGMM));
            
            % [row, col, gmm ] = fem1d_gmm( p, N, P, pot, superM);
            % GMM = sparse( row, col, gmm, N*p+1, N*p+1, length(gmm));

            phi = @(x) x;
            pot = (phi(x1));
            VL = pot(1);
            VR = pot(end);
            LGPMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGPMM = transpose(LGPMM) - diag(diag(LGPMM));

            %% stiffness matrix
            [row1 col1 gsm nnz1_] = fem1d_crgsm( p, N);
            GSM = sparse(row1, col1, gsm, N*p+1, N*p+1, nnz1_);
            %%
            GPMM = -J*(UGPMM+LGPMM); 
            GM_ = GSM/J+(-I*r*J)*(UGMM+LGMM);
            % GM_ = GSM/J+(-I*r*J)*(GMM);
            %%
            u0  = profile(x,0);
            U   = mxfem1d_FLT( N, FM, u0);
            X   = zeros(p*N+1,1);
            Y   = zeros((p+1)*N,1);
            E = U;
            eABC1   =   zeros(1,Nt);
            %%
            YL          =   0;
            phiL        =   0;
            VPhi.p1o2L  =   zeros( 1, 1);
            VPhi.m1o2L  =   zeros( 1, 1);
            VPhi.m1L    =   zeros( 1, 2);

            YR          =   0;
            phiR        =   0;
            VPhi.p1o2R  =   zeros( 1, 1);
            VPhi.m1o2R  =   zeros( 1, 1);
            VPhi.m1R    =   zeros(1,2);

            GM_(1,1)     =   GM_(1,1)+r1o2*ph;
            GM_(N+1,N+1) =   GM_(N+1,N+1)+r1o2*ph;
            spy(GM_);
            %%
            KL2 = KL*KL;
            KR2 = KR*KR;

            EKL2t   =   exp(I*KL2*0.5*dt);
            EKR2t   =   exp(I*KR2*0.5*dt);
            tmpEKL2t   =   exp(-I*KL2*dt);
            tmpEKR2t   =   exp(-I*KR2*dt);

            ir = 1/r;
            Twoir = 2/r;
            tjp1o2 = -0.5*dt;

            phiL0 = chi*AL*conj(AL);
            phiR0 = chi*AR*conj(AR);

            % phr1o2 = ph*r1o2;

            KLAL = KL*AL;
            KRAR = KR*AR;

            phr1o2AL = ph*r1o2*AL;
            phr1o2AR = ph*r1o2*AR;

            iphr1o2pKL = (iph*r1o2+KL);
            iphr1o2mKR = (iph*r1o2-KR);

            vG1o2   =  [];
            vGm1o2  =  [];
            for j=2:Nt,
                Pvb =   mxfem1d_PrjL2F( p, E);
                F       =   -I*r*J*Pvb;
                vG1o2   =   [G1o2(j),vG1o2];
                vGm1o2  =   [Gm1o2(j),vGm1o2];
                % Left:  u_x-i^(-1/2)*r1o2*u=BL
                % right: u_x+i^(-1/2)*r1o2*u=BR
                
                SL(1)   =   ph*sum(vG1o2.*VPhi.p1o2L);
                SL(2)   =   -iph*sum(vGm1o2.*VPhi.m1o2L);
                SR(1)   =   -ph*sum(vG1o2.*VPhi.p1o2R);
                SR(2)   =   iph*sum(vGm1o2.*VPhi.m1o2R);
                Res     =   1;
                it      =   0;
                
                tjp1o2  =   tjp1o2+dt;
                EKL2t   =   EKL2t*tmpEKL2t;
                EKR2t   =   EKR2t*tmpEKR2t;
                EXiL    =   exp(I*XiL(tjp1o2));
                EXiR    =   exp(I*XiR(tjp1o2));
                tmpV    =   V(tjp1o2);
                tmpdXi  =   dXi(tjp1o2);
                GM      =   GPMM*tmpV+GM_;
                
                EXiLEKL2t = EXiL*EKL2t;
                EXiREKR2t = EXiR*EKR2t;
                
                while (Res>tol && it<Nth)
                    PsiL    =   (YL-AL*EXiLEKL2t);
                    qL      =   (phiL0-phiL)*AL*EXiL;
                    qL1     =   VPhi.m1L(1)+ir*qL;
                    
                    PsiR    =   (YR-AR*EXiREKR2t);
                    qR      =   (phiR0-phiR)*AR*EXiR;
                    qR1     =   VPhi.m1R(1)+ir*qR;

                    BL      =   SL(1)+0.5*(phiL+VL*tmpV)*SL(2)...
                                -0.5*iph*(phiL+VL*tmpV)*rm1o2*PsiL...
                                +(iphr1o2pKL)*qL1*EKL2t;
                    BR      =   SR(1)+0.5*(phiR+VR*tmpV)*SR(2)...
                                +0.5*iph*(phiR+VR*tmpV)*rm1o2*PsiR...
                                -(iphr1o2mKR)*qR1*EKR2t;
                    
                    mBL     =   BL+(I*KLAL+I*AL*tmpdXi-phr1o2AL)*EXiLEKL2t;
                    mBR     =   BR+(I*KRAR+I*AR*tmpdXi+phr1o2AR)*EXiREKR2t;
                    
                    tmpu    =   mxfem1d_ILT( N, IMi, Y);
                    tmpNL   =   chi*(tmpu.^2).*conj(tmpu);
                    tmp     =   mxfem1d_FLT(N, FMi, tmpNL);
                    
                    Pvb     =   mxfem1d_PrjL2F( p, tmp);
                    tmpF    =   Pvb*J+F;
                    tmpF(1)     =   tmpF(1)-mBL;
                    tmpF(N+1)   =   tmpF(N+1)+mBR;
                    
                    X       =   GM\tmpF;
                    
                    tmpX    =   mxfem1d_F2L( p, X);
                    Res     =   sqrt(J)*mxfem1d_Norm2( p, N, Y-tmpX);
                    it = it+1;
                    if Res>=tol && it==Nth,
                        error (['Threshold iteration reached at (j=' num2str(j) ') with Res =' num2str(Res) '.']);
                    end
                    Y       =   tmpX;
                    YL      =   IM(1,:)*Y(1:p+1);
                    YR      =   IM(p+1,:)*Y((p+1)*(N-1)+1:(p+1)*(N-1)+p+1);
                    phiL    =   chi*YL*conj(YL);
                    phiR    =   chi*YR*conj(YR);
                end
                E  =   2*Y-E;
               
                PsiL        =   (YL-AL*EXiLEKL2t);
                qL          =   (phiL0-phiL)*AL*EXiL;
                qL1         =   VPhi.m1L(1)+ir*qL;
                
                VPhi.m1L(1) =   VPhi.m1L(1)+Twoir*qL;
                VPhi.p1o2L(1,j)  =   (PsiL+I*qL1*EKL2t);
                VPhi.m1o2L(1,j)  =   PsiL;
                
                PsiR        =   (YR-AR*EXiREKR2t);
                qR          =   (phiR0-phiR)*AR*EXiR;
                qR1         =   VPhi.m1R(1)+ir*qR;
                
                VPhi.m1R(1) =   VPhi.m1R(1)+Twoir*qR;
                VPhi.p1o2R(1,j)  =   (PsiR+I*qR1*EKR2t);
                VPhi.m1o2R(1,j)  =   PsiR;
                
                u           =   profile(x,(j-1)*dt);
                U           =   mxfem1d_FLT( N, FM, u);
                norm_ana    =   mxfem1d_Norm2( p, N, U);
                eABC1(j)    =   mxfem1d_Norm2( p, N, E-U)/norm_ana;
                disp(eABC1(j))
            end
            out = eABC1;
        end
        
        function out = ABC1a_Pade1(Data)
            % Solving GP equation using FEM
            % iu_t+u_xx+xg(t)u+chi*|u|^2u=0
            % with artificial boundary conditions ABC_1.
            % ABC_1^([M/M]) -- ABC_1 discretized using Pade approximation.
            % Argument structure is 
            % Data
            % |----fp  --->  p,     Highest degree of the polynomials
            % |----fN  --->  N,     Number of partitions defining finite elements
            % |----fx  --->  x,     Grid points
            % |----fFM --->  FM,    Forward Transform matrix
            % |----fIM --->  IM,    Inverse transform matrix
            % |----fFMi--->  FMi,   Forward transform for anti-aliasing
            % |----fIMi--->  IMi,   Inverse transform for anti-aliasing
            % |----fP  --->  P,     Number of LGL points for Numerical integration
            % |----fxri--->  xri,   LGL points for Numerical integration
            % |----fquadWi---> quadWi,    Quadrature weights for LGL points
            % |----ftol--->  tol,   tolerance
            % |----fdt --->  dt,    Time step
            % |----fNt --->  Nt,    Number of time steps
            % |----fNth--->  Nth,   Threshold iteration
            % |----fMp --->  Mp,    Order of the Pade approximants
            % |----fx_r--->  x_r,   Left boundary
            % |----fx_l--->  x_l,   Right boundary
            % |----fAL --->  AL,    Value of the field on the left exterior domain
            % |----fAR --->  AR,    Value of the field on the right exterior domain
            % |----fKL --->  KL,    Wave vector for the left exterior domain
            % |----fKR --->  KR,    Wave vector for the right exterior domain
            % |----fchi--->  chi,   Nonlinear coefficient (refer to the evolution equation)
            % |----fu0 --->  u0,    Initial data
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %
            % Degree of the Legendre polynomials involved --> p
            % tol = 1e-9;
            % [xr quadW IM FM] = LGLdataLT2( p, tol);
            % xr: LGL points
            % FM: matrix for forward transform
            % IM: matrix for inverse transform

            % Anti aliasing for the nonlinear term
            % P = 2*p;% For the purpose of numerical integration
            % Degree of the Legendre polynomials involved --> P-1
            % [xri quadWi] = LGLdataLT1( P, tol);
            % FMi = LGLdataFM( p, P, xri); 
            % IMi = LGLdataIM( P, p, xri);
            %
            % Pade method for fractional derivativesm requires 
            % the mex function PF_Pade.
            % 
            % Date : Apr 2012
            %%
            FM      =   Data.fFM;
            IM      =   Data.fIM;
            FMi     =   Data.fFMi;
            IMi     =   Data.fIMi;
            quadWi  =   Data.fquadWi;
            xri     =   Data.fxri;
            tol     =   Data.ftol;
            N       =   Data.fN;
            p       =   Data.fp;
            P       =   Data.fP;
            dt      =   Data.fdt;
            Nt      =   Data.fNt;
            Nth     =   Data.fNth;
            Mp      =   Data.fMp;
            x_r     =   Data.fx_r;
            x_l     =   Data.fx_l;
            AL      =   Data.fAL;
            AR      =   Data.fAR;
            KL      =   Data.fKL;
            KR      =   Data.fKR;
            chi     =   Data.fchi;
            x       =   Data.fx;
            V       =   str2func(Data.hpot);
            XiL     =   str2func(Data.hXiL);
            XiR     =   str2func(Data.hXiR);
            dXi     =   str2func(Data.hdXi);
            profile =   str2func(Data.hfunc);


            I   =   sqrt(-1);
            dx  = (x_r-x_l)/N;
            J   = dx/2;
            r   =   2/dt;
             
            ph  =   exp(-I*pi/4);
            iph =   exp(I*pi/4);
            %% Pade approximation
            nu  = 1/2;
            [Num Denom] = pf_pade_exp( Mp, nu, tol);
            b0          =   Num(1);
            b1o2        =   Num(2:end);
            G1o2        =   r./(r+Denom);
            r1o2        =   (b0+sum(b1o2.*G1o2)/r);
            bG1o2       =   2*b1o2.*G1o2/r;

            nu  = -1/2;
            [Num Denom] = pf_pade_exp( Mp, nu, tol);
            b0          =   Num(1);
            bm1o2        =   Num(2:end);
            Gm1o2        =   r./(r+Denom);
            rm1o2        =   (b0+sum(bm1o2.*Gm1o2)/r);
            bGm1o2       =   2*bm1o2.*Gm1o2/r;
            %%
            superM = fem1d_superM( p, P, IMi, quadWi);
            %phi = @(x) 1+x.^2;
            x1 = ref2mesh (xri, N, [x_l x_r]);
            % mass matrix
            pot = ones(size(x1));
            [row col nnz_] = fem1d_crgmmind( p, N);
            [ugmm lgmm] = fem1d_gmm( p, N, P, pot, superM, nnz_);
            UGMM = sparse( row, col, ugmm, N*p+1, N*p+1, nnz_);
            LGMM = sparse( col, row, lgmm, N*p+1, N*p+1, nnz_);
            phi = @(x) x;
            pot = (phi(x1));
            VL = pot(1);
            VR = pot(end);
            [ugpmm lgpmm] = fem1d_gmm( p, N, P, pot, superM, nnz_);
            UGPMM = sparse( row, col, ugpmm, N*p+1, N*p+1, nnz_);
            LGPMM = sparse( col, row, lgpmm, N*p+1, N*p+1, nnz_);

            %% stiffness matrix
            [row1 col1 gsm nnz1_] = fem1d_crgsm( p, N);
            GSM = sparse(row1, col1, gsm, N*p+1, N*p+1, nnz1_);
            %%
            GPMM = -J*(UGPMM+LGPMM); 
            GM_ = GSM/J+(-I*r*J)*(UGMM+LGMM);
            %%
            u0  =   profile(x,0);
            U   = fem_flt( p, N, p, FM, u0);
            X   = zeros(p*N+1,1);
            Y   = zeros((p+1)*N,1);
            E = U;
            eABC1   =   zeros(1,Nt);
            %%
            YL          =   0;
            phiL        =   0;
            VPhi.p1o2L  =   zeros( 1, Mp);
            VPhi.m1o2L  =   zeros( 1, Mp);
            VPhi.m1L    =   zeros( 1, 2);

            YR          =   0;
            phiR        =   0;
            VPhi.p1o2R  =   zeros( 1, Mp);
            VPhi.m1o2R  =   zeros( 1, Mp);
            VPhi.m1R    =   zeros(1,2);

            GM_(1,1)     =   GM_(1,1)+r1o2*ph;
            GM_(N+1,N+1) =   GM_(N+1,N+1)+r1o2*ph;
            % spy(GM);
            %%
            KL2 = KL*KL;
            KR2 = KR*KR;

            EKL2t   =   exp(I*KL2*0.5*dt);
            EKR2t   =   exp(I*KR2*0.5*dt);
            tmpEKL2t   =   exp(-I*KL2*dt);
            tmpEKR2t   =   exp(-I*KR2*dt);

            ir = 1/r;
            Twoir = 2/r;
            tjp1o2 = -0.5*dt;

            phiL0 = chi*AL*conj(AL);
            phiR0 = chi*AR*conj(AR);

            % phr1o2 = ph*r1o2;

            KLAL = KL*AL;
            KRAR = KR*AR;

            phr1o2AL = ph*r1o2*AL;
            phr1o2AR = ph*r1o2*AR;

            iphr1o2pKL = (iph*r1o2+KL);
            iphr1o2mKR = (iph*r1o2-KR);

            TwoG1o2m1 = (2*G1o2-1);
            TwoGm1o2m1 = (2*Gm1o2-1);

            for j=2:Nt,
                [Pv Pb] =   prjlp2fem( p, N, E);
                F       =   -I*r*J*[Pv;Pb];
                % Left:  u_x-i^(-1/2)*r1o2*u=BL
                % right: u_x+i^(-1/2)*r1o2*u=BR
                
                SL(1)   =   ph*sum(G1o2.*VPhi.p1o2L);
                SL(2)   =   -iph*sum(Gm1o2.*VPhi.m1o2L);
                SR(1)   =   -ph*sum(G1o2.*VPhi.p1o2R);
                SR(2)   =   iph*sum(Gm1o2.*VPhi.m1o2R);
                Res     =   1;
                it      =   0;
                
                tjp1o2  =   tjp1o2+dt;
                EKL2t   =   EKL2t*tmpEKL2t;
                EKR2t   =   EKR2t*tmpEKR2t;
                EXiL    =   exp(I*XiL(tjp1o2));
                EXiR    =   exp(I*XiR(tjp1o2));
                tmpV    =   V(tjp1o2);
                tmpdXi  =   dXi(tjp1o2);
                GM      =   GPMM*tmpV+GM_;
                
                EXiLEKL2t = EXiL*EKL2t;
                EXiREKR2t = EXiR*EKR2t;
                
                while (Res>tol && it<Nth)
                    PsiL    =   (YL-AL*EXiLEKL2t);
                    qL      =   (phiL0-phiL)*AL*EXiL;
                    qL1     =   VPhi.m1L(1)+ir*qL;
                    
                    PsiR    =   (YR-AR*EXiREKR2t);
                    qR      =   (phiR0-phiR)*AR*EXiR;
                    qR1     =   VPhi.m1R(1)+ir*qR;

                    BL      =   SL(1)+0.5*(phiL+VL*tmpV)*SL(2)...
                                -0.5*iph*(phiL+VL*tmpV)*rm1o2*PsiL...
                                +(iphr1o2pKL)*qL1*EKL2t;
                    BR      =   SR(1)+0.5*(phiR+VR*tmpV)*SR(2)...
                                +0.5*iph*(phiR+VR*tmpV)*rm1o2*PsiR...
                                -(iphr1o2mKR)*qR1*EKR2t;
                    
                    mBL     =   BL+(I*KLAL+I*AL*tmpdXi-phr1o2AL)*EXiLEKL2t;
                    mBR     =   BR+(I*KRAR+I*AR*tmpdXi+phr1o2AR)*EXiREKR2t;
                    
                    tmpu    =   fem_ilt( p, N, P, IMi, Y);
                    tmpNL   =   chi*(tmpu.^2).*conj(tmpu);
                    tmp     =   fem_flt(p, N, P, FMi, tmpNL);
                    
                    [Pv Pb] =   prjlp2fem( p, N, tmp);
                    tmpF    =   [Pv;Pb]*J+F;
                    tmpF(1)     =   tmpF(1)-mBL;
                    tmpF(N+1)   =   tmpF(N+1)+mBR;
                    
                    X       =   GM\tmpF;
                    
                    tmpX    =   fem2lp( p, N, X(1:N+1), X(N+2:end));
                    Res     =   sqrt(J)*femlp_norm2( p, N, Y-tmpX);
                    it = it+1;
                    if Res>=tol && it==Nth,
                        error (['Threshold iteration reached at (j=' num2str(j) ') with Res =' num2str(Res) '.']);
                    end
                    Y       =   tmpX;
                    YL      =   IM(1,:)*Y(1:p+1);
                    YR      =   IM(p+1,:)*Y((p+1)*(N-1)+1:(p+1)*(N-1)+p+1);
                    phiL    =   chi*YL*conj(YL);
                    phiR    =   chi*YR*conj(YR);
                end
                E  =   2*Y-E;
               
                PsiL        =   (YL-AL*EXiLEKL2t);
                qL          =   (phiL0-phiL)*AL*EXiL;
                qL1         =   VPhi.m1L(1)+ir*qL;
                
                VPhi.m1L(1) =   VPhi.m1L(1)+Twoir*qL;
                VPhi.p1o2L  =   TwoG1o2m1.*VPhi.p1o2L+bG1o2*(PsiL+I*qL1*EKL2t);
                VPhi.m1o2L  =   TwoGm1o2m1.*VPhi.m1o2L+bGm1o2*PsiL;
                
                PsiR        =   (YR-AR*EXiREKR2t);
                qR          =   (phiR0-phiR)*AR*EXiR;
                qR1         =   VPhi.m1R(1)+ir*qR;
                
                VPhi.m1R(1) =   VPhi.m1R(1)+Twoir*qR;
                VPhi.p1o2R  =   TwoG1o2m1.*VPhi.p1o2R+bG1o2*(PsiR+I*qR1*EKR2t);
                VPhi.m1o2R  =   TwoGm1o2m1.*VPhi.m1o2R+bGm1o2*PsiR;
                
                u           =   profile(x,(j-1)*dt);
                U           =   fem_flt( p, N, p, FM, u);
                norm_ana    =   femlp_norm2( p, N, U);
                eABC1(j)    =   femlp_norm2( p, N, E-U)/norm_ana;
            end
            out = eABC1;
        end

        function r = gpe_bs( xx, tt, a0, eta, c, g_0,g_1, mu)
            [x,t] = ndgrid(xx,tt);
            beta  = g_0*t + (g_1*sin(mu*t))/mu;
            nu    = (2*g_1*cos(mu*t))/mu^2 - g_0*t.^2 - (2*g_1)/mu^2;

            ibeta2 = (mu*((g_1^2*t)/2 - 2*g_0*g_1*t.*cos(mu*t)) -... 
                (g_1^2*sin(2*mu*t))/4 + 2*g_0*g_1*sin(mu*t))/mu^3 + (g_0^2*t.^3)/3;

            Xi = (c/2)*nu+beta.*x-ibeta2+2*(a0^2)*t;

            eta1 = sqrt(eta*eta+4*a0*a0);
            p    = eta*eta1;
            E    = exp(1i*(c*x/2-c^2*t/4+Xi));

            un = eta*(eta*cos(p*t)+1i*eta1*sin(p*t));
            ud = (eta1*cosh(eta*(x+nu-c*t))+2*a0*cos(p*t));
            r  = (a0+un./ud).*E;
        end
        
        function r = Xi(x,t,a0,k,g_0,g_1,mu,chi)
            c = 2*k;
            beta = g_0*t + (g_1*sin(mu*t))/mu;
            nu = (2*g_1*cos(mu*t))/mu^2 - g_0*t.^2 - (2*g_1)/mu^2;

            ibeta2 = (mu*((g_1^2*t)/2 - 2*g_0*g_1*t*cos(mu*t)) -... 
                (g_1^2*sin(2*mu*t))/4 + 2*g_0*g_1*sin(mu*t))/mu^3 + (g_0^2*t.^3)/3;
             
            r = (c/2)*nu+beta.*x-ibeta2+chi*(a0)*conj(a0)*t;
        end

        function beta = dXi(t,g_0,g_1,mu)
            beta = g_0*t+g_1*sin(mu*t)/mu;
        end
    end
end
