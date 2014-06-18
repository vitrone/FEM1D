classdef LSE
    
    properties(Constant)
        CQ_M0 = 20;
    
    end


    methods (Static)
        
        function wh = quadW_CQ_BDF1(N, dx, dt, tol )
            if(N<LSE.CQ_M0)
                N = LSE.CQ_M0;
            end
            F  = @(s) exp(-dx*sqrt(-1i*(s))).*sqrt(-1i*s);
            sh = @(zeta, dt) (1/dt)*(1-zeta);
            
            n    = (0:N-1);
            rho  = (tol)^(0.5/N);
            zeta = rho*exp(1i*2*pi*n/N);
            Sh   = sh(zeta, dt);
            
            F_val = F(Sh);
            % Quadrature weights
            wh = (rho.^(-n)).*fft(F_val)/N;
        end

        function wh = quadW_CQ_BDF2(N, dx, dt, tol )
            if(N<LSE.CQ_M0)
                N = LSE.CQ_M0;
            end
            F  = @(s) exp(-dx*sqrt(-1i*(s))).*sqrt(-1i*s);
            sh = @(zeta, dt) (1/dt)*(3/2-2*zeta+zeta.^2/2);
            
            n    = (0:N-1);
            rho  = (tol)^(0.5/N);
            zeta = rho*exp(1i*2*pi*n/N);
            Sh   = sh(zeta, dt);
            
            F_val = F(Sh);
            % Quadrature weights
            wh = (rho.^(-n)).*fft(F_val)/N;
        end

        % Trapezoidal rule
        function wh = quadW_CQ(N, dx, dt, tol )
            if(N<LSE.CQ_M0)
                N = LSE.CQ_M0;
            end
            F  = @(s) exp(-dx*sqrt(-1i*(s))).*sqrt(-1i*s);
            sh = @(zeta, dt) (2/dt)*(1-zeta)./(1+zeta);
            
            n    = (0:N-1);
            rho  = (tol)^(0.5/N);
            zeta = rho*exp(1i*2*pi*n/N);
            Sh   = sh(zeta, dt);
            
            F_val = F(Sh);
            % Quadrature weights
            wh = (rho.^(-n)).*fft(F_val)/N;
        end

        function wh = quadW_CQ2(N, dx, dt, tol )
            if(N<LSE.CQ_M0)
                N = LSE.CQ_M0;
            end
            F  = @(s) exp(-dx*sqrt(-1i*(s)))./sqrt(-1i*s);
            sh = @(zeta, dt) (2/dt)*(1-zeta)./(1+zeta);
            
            n    = (0:N-1);
            rho  = (tol)^(0.5/N);
            zeta = rho*exp(1i*2*pi*n/N);
            Sh   = sh(zeta, dt);
            
            F_val = F(Sh);
            % Quadrature weights
            wh = (rho.^(-n)).*fft(F_val)/N;
        end


        function E = TBC_CQ_evol(Data)
            % Solving Linear Schroedinger equation using FEM
            % iu_t+u_xx=0
            % TBC^(*) -- discretized using convolution quadrature (CQ).
            % Argument structure is 
            % Data
            %   |----p  ---->   Highest degree of the polynomials
            %   |----N  ---->   Number of partitions defining finite elements
            %   |----x  ---->   Grid points
            %   |----FM ---->   Forward Transform matrix
            %   |----IM ---->   Inverse transform matrix
            %   |----FMi---->   Forward transform for anti-aliasing
            %   |----IMi---->   Inverse transform for anti-aliasing
            %   |----P  ---->   Number of LGL points for Numerical integration
            %   |----xri---->   LGL points for Numerical integration
            %   |----quadWi->   Quadrature weights for LGL points
            %   |----tol---->   tolerance
            %   |----dt ---->   Time step
            %   |----Nt ---->   Number of time steps
            %   |----Nth---->   Threshold iteration
            %   |----x_r---->   Left boundary
            %   |----x_l---->   Right boundary
            %   |----u0 ---->   Initial data
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %
            FM      = Data.FM;
            IM      = Data.IM;
            FMi     = Data.FMi;
            IMi     = Data.IMi;
            quadWi  = Data.quadWi;
            xri     = Data.xri;
            tol     = Data.tol;
            N       = Data.N;
            p       = Data.p;
            P       = Data.P;
            dt      = Data.dt;
            Nt      = Data.Nt;
            x_r     = Data.x_r;
            x_l     = Data.x_l;
            dx_l    = Data.dx_l;
            dx_r    = Data.dx_r;
            x       = Data.x;
            tol     = Data.tol;
            profile = Data.hfunc;

            I  =  sqrt(-1);
            dx =  (x_r-x_l)/N;
            J  =  dx/2;
            r  =  2/dt;

            Q = mxfem1d_quadM( quadWi, IMi);
            x1 = ref2mesh (xri, N, [x_l x_r]);

            % mass matrix
            pot = ones(size(x1));
            LGMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGMM = transpose(LGMM) - diag(diag(LGMM));
            
            %% stiffness matrix
            [row1, col1, gsm, nnz1_] = fem1d_crgsm( p, N);
            GSM = sparse(row1, col1, gsm, N*p+1, N*p+1, nnz1_);
            %%
            GM = GSM/J+(-I*r*J)*(UGMM+LGMM);
            u0 = profile(x,0);
            U  = mxfem1d_FLT( N, FM, u0);
            X  = zeros(p*N+1,1);
            Y  = zeros((p+1)*N,1);

            E = zeros(N*(p+1), Nt);
            E(:, 1) = U;

            yL = [0];
            yR = [0];

            SL = 0;
            SR = 0;

            whr = LSE.quadW_CQ(Nt, dx_r, dt, tol);
            whl = LSE.quadW_CQ(Nt, dx_l, dt, tol);

            if(dx_r==0)
                nx_l  = 0; 
                xi_l  = -1;
                IMb_l = mxLGLdataIM( p, xi_l);
                GM(1, 1) = GM(1, N)+whl(1);
                nx_r  = N-1; 
                xi_r  = 1;
                IMb_r = mxLGLdataIM( p, xi_r);
                GM(N+1, N+1) = GM(N+1, N+1)+whr(1);
            else
                [nx_l, xi_l] = LSE.fem_index(N, x_l, x_r, x_l+dx_l);
                IMb_l   = mxLGLdataIM( p, xi_l);
                IMvb_l  = LSE.Lobatto(IMb_l);
                GM(1, nx_l+1) = GM(1, nx_l+1) + whl(1)*IMvb_l(1);
                GM(1, nx_l+2) = GM(1, nx_l+2) + whl(1)*IMvb_l(2);
                shift = N+1 + nx_l*(p-2);
                for k=3:p+1
                    GM(1, shift+k) = GM(1, shift+k) + ...
                                       whl(1)*IMvb_l(k);
                end

                [nx_r, xi_r] = LSE.fem_index(N, x_l, x_r, x_r-dx_r);
                IMb_r  = mxLGLdataIM( p, xi_r);
                IMvb_r = LSE.Lobatto(IMb_r);
                GM(N+1, nx_r+1) = GM(N+1, nx_r+1) + whr(1)*IMvb_r(1);
                GM(N+1, nx_r+2) = GM(N+1, nx_r+2) + whr(1)*IMvb_r(2);
                shift = N+1 + nx_r*(p-2);
                for k=3:p+1
                    GM(N+1, shift+k) = GM(N+1, shift+k) + ...
                                       whr(1)*IMvb_r(k);
                end
            end
            
            for j=2:Nt,
                Pvb = mxfem1d_PrjL2F( p, E(:, j-1));
                F   = -I*r*J*Pvb;

                whl_tmp = fliplr(whl((2:j)));
                whr_tmp = fliplr(whr((2:j)));

                SL = sum(whl_tmp.*yL);
                SR = sum(whr_tmp.*yR);

                tmpF = F;
                tmpF(1) = F(1)-SL;
                tmpF(N+1) = F(N+1)-SR;

                X = GM\tmpF;
                Y = mxfem1d_F2L( p, X);

                E(:, j)  = 2*Y-E(:, j-1);
                yL(1, j) = IMb_l*Y((nx_l)*(p+1)+1:(nx_l)*(p+1)+p+1,1);
                yR(1, j) = IMb_r*Y((nx_r)*(p+1)+1:(nx_r)*(p+1)+p+1,1);
            end
        end

        function eTBC = TBC_CQ_error(Data)
            % Solving Linear Schroedinger equation using FEM
            % iu_t+u_xx=0
            % TBC^(*) -- discretized using convolution quadrature (CQ).
            % Argument structure is 
            % Data
            %   |----p  ---->   Highest degree of the polynomials
            %   |----N  ---->   Number of partitions defining finite elements
            %   |----x  ---->   Grid points
            %   |----FM ---->   Forward Transform matrix
            %   |----IM ---->   Inverse transform matrix
            %   |----FMi---->   Forward transform for anti-aliasing
            %   |----IMi---->   Inverse transform for anti-aliasing
            %   |----P  ---->   Number of LGL points for Numerical integration
            %   |----xri---->   LGL points for Numerical integration
            %   |----quadWi->   Quadrature weights for LGL points
            %   |----tol---->   tolerance
            %   |----dt ---->   Time step
            %   |----Nt ---->   Number of time steps
            %   |----Nth---->   Threshold iteration
            %   |----x_r---->   Left boundary
            %   |----x_l---->   Right boundary
            %   |----u0 ---->   Initial data
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %
            FM      = Data.FM;
            IM      = Data.IM;
            FMi     = Data.FMi;
            IMi     = Data.IMi;
            quadWi  = Data.quadWi;
            xri     = Data.xri;
            tol     = Data.tol;
            N       = Data.N;
            p       = Data.p;
            P       = Data.P;
            dt      = Data.dt;
            Nt      = Data.Nt;
            x_r     = Data.x_r;
            x_l     = Data.x_l;
            dx_l    = Data.dx_l;
            dx_r    = Data.dx_r;
            x       = Data.x;
            tol     = Data.tol;
            profile = Data.hfunc;

            I  =  sqrt(-1);
            dx =  (x_r-x_l)/N;
            J  =  dx/2;
            r  =  2/dt;

            Q = mxfem1d_quadM( quadWi, IMi);
            x1 = ref2mesh (xri, N, [x_l x_r]);

            % mass matrix
            pot = ones(size(x1));
            LGMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGMM = transpose(LGMM) - diag(diag(LGMM));
            
            %% stiffness matrix
            [row1, col1, gsm, nnz1_] = fem1d_crgsm( p, N);
            GSM = sparse(row1, col1, gsm, N*p+1, N*p+1, nnz1_);
            %%
            GM = GSM/J+(-I*r*J)*(UGMM+LGMM);
            u0 = profile(x,0);
            U  = mxfem1d_FLT( N, FM, u0);
            X  = zeros(p*N+1,1);
            Y  = zeros((p+1)*N,1);

            E = zeros(N*(p+1), Nt);
            E = U;

            yL = [0];
            yR = [0];

            SL = 0;
            SR = 0;

            whr = LSE.quadW_CQ(Nt, dx_r, dt, tol);
            whl = LSE.quadW_CQ(Nt, dx_l, dt, tol);

            if(dx_r==0)
                nx_l  = 0; 
                xi_l  = -1;
                IMb_l = mxLGLdataIM( p, xi_l);
                GM(1, 1) = GM(1, 1)+whl(1);
                nx_r  = N-1; 
                xi_r  = 1;
                IMb_r = mxLGLdataIM( p, xi_r);
                GM(N+1, N+1) = GM(N+1, N+1)+whr(1);
            else
                [nx_l, xi_l] = LSE.fem_index(N, x_l, x_r, x_l+dx_l);
                IMb_l   = mxLGLdataIM( p, xi_l);
                IMvb_l  = LSE.Lobatto(IMb_l);
                GM(1, nx_l+1) = GM(1, nx_l+1) + whl(1)*IMvb_l(1);
                GM(1, nx_l+2) = GM(1, nx_l+2) + whl(1)*IMvb_l(2);
                shift = N+1 + nx_l*(p-2);
                for k=3:p+1
                    GM(1, shift+k) = GM(1, shift+k) + ...
                                       whl(1)*IMvb_l(k);
                end

                [nx_r, xi_r] = LSE.fem_index(N, x_l, x_r, x_r-dx_r);
                IMb_r  = mxLGLdataIM( p, xi_r);
                IMvb_r = LSE.Lobatto(IMb_r);
                GM(N+1, nx_r+1) = GM(N+1, nx_r+1) + whr(1)*IMvb_r(1);
                GM(N+1, nx_r+2) = GM(N+1, nx_r+2) + whr(1)*IMvb_r(2);
                shift = N+1 + nx_r*(p-2);
                for k=3:p+1
                    GM(N+1, shift+k) = GM(N+1, shift+k) + ...
                                       whr(1)*IMvb_r(k);
                end
            end
            
            for j=2:Nt,
                Pvb = mxfem1d_PrjL2F( p, E);
                F   = -I*r*J*Pvb;

                whl_tmp = fliplr(whl((2:j)));
                whr_tmp = fliplr(whr((2:j)));

                SL = sum(whl_tmp.*yL);
                SR = sum(whr_tmp.*yR);

                tmpF = F;
                tmpF(1) = F(1)-SL;
                tmpF(N+1) = F(N+1)-SR;

                X = GM\tmpF;
                Y = mxfem1d_F2L( p, X);

                E = 2*Y-E;
                yL(1, j) = IMb_l*Y((nx_l)*(p+1)+1:(nx_l)*(p+1)+p+1,1);
                yR(1, j) = IMb_r*Y((nx_r)*(p+1)+1:(nx_r)*(p+1)+p+1,1);

                u = profile(x,(j-1)*dt);
                U = mxfem1d_FLT( N, FM, u);
                norm_ana = mxfem1d_Norm2( p, N, U);
                eTBC(j)  = mxfem1d_Norm2( p, N, E-U)/norm_ana;
            end
            out = eTBC;
        end

        function [eTBC, E] = TBC_CQ_evol3(Data)
            % Solving Linear Schroedinger equation using FEM
            % iu_t+u_xx=0
            % TBC^(*) -- discretized using convolution quadrature (CQ).
            % Argument structure is 
            % Data
            %   |----p  ---->   Highest degree of the polynomials
            %   |----N  ---->   Number of partitions defining finite elements
            %   |----x  ---->   Grid points
            %   |----FM ---->   Forward Transform matrix
            %   |----IM ---->   Inverse transform matrix
            %   |----FMi---->   Forward transform for anti-aliasing
            %   |----IMi---->   Inverse transform for anti-aliasing
            %   |----P  ---->   Number of LGL points for Numerical integration
            %   |----xri---->   LGL points for Numerical integration
            %   |----quadWi->   Quadrature weights for LGL points
            %   |----tol---->   tolerance
            %   |----dt ---->   Time step
            %   |----Nt ---->   Number of time steps
            %   |----Nth---->   Threshold iteration
            %   |----x_r---->   Left boundary
            %   |----x_l---->   Right boundary
            %   |----u0 ---->   Initial data
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %
            FM      = Data.FM;
            IM      = Data.IM;
            FMi     = Data.FMi;
            IMi     = Data.IMi;
            quadWi  = Data.quadWi;
            xri     = Data.xri;
            tol     = Data.tol;
            N       = Data.N;
            p       = Data.p;
            P       = Data.P;
            dt      = Data.dt;
            Nt      = Data.Nt;
            x_r     = Data.x_r;
            x_l     = Data.x_l;
            dx_l    = Data.dx_l;
            dx_r    = Data.dx_r;
            x       = Data.x;
            profile = (Data.hfunc);
            V       = (Data.V);
            nu      = (Data.nu);
            beta    = (Data.beta);
            B       = (Data.B);

            I  =  sqrt(-1);
            dx =  (x_r-x_l)/N;
            J  =  dx/2;
            r  =  2/dt;

            Q = mxfem1d_quadM( quadWi, IMi);
            x1 = ref2mesh (xri, N, [x_l x_r]);

            % mass matrix
            pot = ones(size(x1));
            LGMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGMM = transpose(LGMM) - diag(diag(LGMM));
            
            %% stiffness matrix
            [row1, col1, gsm, nnz1_] = fem1d_crgsm( p, N);
            GSM = sparse(row1, col1, gsm, N*p+1, N*p+1, nnz1_);

            phi = @(x) x;
            pot = (phi(x1));
            LGPMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGPMM = transpose(LGPMM) - diag(diag(LGPMM));
            GPMM  = -J*(UGPMM+LGPMM);  
            %%
            GM = GSM/J+(-I*r*J)*(UGMM+LGMM);
            u0 = profile(x,0);
            U  = mxfem1d_FLT( N, FM, u0);
            X  = zeros(p*N+1,1);
            Y  = zeros((p+1)*N,1);

            E = zeros(N*(p+1), Nt);
            E(:, 1) = U;

            eTBC = zeros(1,Nt);
            tjp1o2 = -0.5*dt;
            for j=2:Nt,

                tjp1o2 = tjp1o2 + dt;

                tmpV = V(tjp1o2);
                GM1 = GPMM*tmpV + GM;

                Pvb = mxfem1d_PrjL2F( p, E(:, j-1));
                F   = -I*r*J*Pvb;

                X = GM1\F;
                Y = mxfem1d_F2L( p, X);

                E(:, j)  = 2*Y-E(:, j-1);


                u = profile(x,(j-1)*dt);
                U = mxfem1d_FLT( N, FM, u);
                norm_ana = mxfem1d_Norm2( p, N, U);
                eTBC(j)  = mxfem1d_Norm2( p, N, E(:, j)-U)/norm_ana;
            end
        end

        function [eTBC, E] = TBC_CQ_evol4(Data)
            % Solving Linear Schroedinger equation using FEM
            % iu_t+u_xx=0
            % TBC^(*) -- discretized using convolution quadrature (CQ).
            % Argument structure is 
            % Data
            %   |----p  ---->   Highest degree of the polynomials
            %   |----N  ---->   Number of partitions defining finite elements
            %   |----x  ---->   Grid points
            %   |----FM ---->   Forward Transform matrix
            %   |----IM ---->   Inverse transform matrix
            %   |----FMi---->   Forward transform for anti-aliasing
            %   |----IMi---->   Inverse transform for anti-aliasing
            %   |----P  ---->   Number of LGL points for Numerical integration
            %   |----xri---->   LGL points for Numerical integration
            %   |----quadWi->   Quadrature weights for LGL points
            %   |----tol---->   tolerance
            %   |----dt ---->   Time step
            %   |----Nt ---->   Number of time steps
            %   |----Nth---->   Threshold iteration
            %   |----x_r---->   Left boundary
            %   |----x_l---->   Right boundary
            %   |----u0 ---->   Initial data
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %
            FM      = Data.FM;
            IM      = Data.IM;
            FMi     = Data.FMi;
            IMi     = Data.IMi;
            quadWi  = Data.quadWi;
            xri     = Data.xri;
            tol     = Data.tol;
            N       = Data.N;
            p       = Data.p;
            P       = Data.P;
            dt      = Data.dt;
            Nt      = Data.Nt;
            x_r     = Data.x_r;
            x_l     = Data.x_l;
            dx_l    = Data.dx_l;
            dx_r    = Data.dx_r;
            x       = Data.x;
            profile = (Data.hfunc);
            V       = (Data.V);
            nu      = (Data.nu);
            beta    = (Data.beta);
            B       = (Data.B);

            I  =  sqrt(-1);
            dx =  (x_r-x_l)/N;
            J  =  dx/2;
            r  =  2/dt;

            Q = mxfem1d_quadM( quadWi, IMi);
            x1 = ref2mesh (xri, N, [x_l x_r]);

            % mass matrix
            pot = ones(size(x1));
            LGMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGMM = transpose(LGMM) - diag(diag(LGMM));
            
            %% stiffness matrix
            [row1, col1, gsm, nnz1_] = fem1d_crgsm( p, N);
            GSM = sparse(row1, col1, gsm, N*p+1, N*p+1, nnz1_);

            phi = @(x) x;
            pot = (phi(x1));
            LGPMM = mxfem1d_sparse_GMM( p, Q, pot);
            UGPMM = transpose(LGPMM) - diag(diag(LGPMM));
            GPMM  = -J*(UGPMM+LGPMM);  
            %%
            GM = GSM/J+(-I*r*J)*(UGMM+LGMM);
            u0 = profile(x,0);
            U  = mxfem1d_FLT( N, FM, u0);
            X  = zeros(p*N+1,1);
            Y  = zeros((p+1)*N,1);

            E = zeros(N*(p+1), Nt);
            E(:, 1) = U;

            yL  = [0];
            yR1 = [0];

            SL = 0;
            SR1 = 0;

            tol = 1e-9;

            tjp1o2 = -0.5*dt;
            t = 0;

            eTBC = zeros(1,Nt);
            wh0_old = 0;
            SR1_old = 0;
            for j=2:Nt,

                t = t + dt;
                nu_tmp1   = nu(t);
                beta_tmp1 = beta(t);
                B_tmp1    = B(t);

                tjp1o2 = tjp1o2 + dt;
                nu_tmp   = nu(tjp1o2);
                beta_tmp = beta(tjp1o2);
                B_tmp    = B(tjp1o2);
                wh1 = LSE.quadW_CQ(j, dx_r + nu_tmp1, dt, tol);

                tmpV = V(tjp1o2);
                GM1 = GPMM*tmpV + GM;

                x_tmp1 = x_r - dx_r - nu_tmp1;

                if(x_tmp1>x_r)
                    err_str = sprintf( 'x-coordinate inconsistent: %0.5f',...
                                       x_tmp1)
                    error(err_str)
                end

                [nx, xi] = LSE.fem_index(N, x_l, x_r, x_tmp1);

                IM_tbc   = mxLGLdataIM( p, xi);
                 
                wh0_new = wh1(1)*exp(1i*beta_tmp1*(dx_r + nu_tmp1));

                IM1_tbc = 0.5*(wh0_new+wh0_old)*LSE.Lobatto(IM_tbc);

                GM1(N+1, N+1)   = GM1(N+1, N+1)  - 1i*beta_tmp;
                GM1(N+1, nx+1)  = GM1(N+1, nx+1) + IM1_tbc(1);
                GM1(N+1, nx+2)  = GM1(N+1, nx+2) + IM1_tbc(2);
                
                shift = N + 1 + nx*(p-2);
                
                for k=3:p+1
                    GM1(N+1, shift+k) = GM1(N+1, shift+k) + ...
                                      IM1_tbc(k);
                end


                Pvb = mxfem1d_PrjL2F( p, E(:, j-1));
                F   = -I*r*J*Pvb;

                wh1_tmp = fliplr(wh1((2:j)));
                Xi1 = beta_tmp1*x_tmp1 - B_tmp1; 
                Xi = beta_tmp1*x_r - B_tmp1; 
                SR1_new = exp(1i*Xi)*sum(wh1_tmp.*yR1);


                tmpF = F;
                tmpF(N+1) = F(N+1)-0.5*(SR1_new + SR1_old);

                X = GM1\tmpF;
                Y = mxfem1d_F2L( p, X);

                E(:, j)  = 2*Y-E(:, j-1);
                yR1(1, j) = exp(-1i*Xi1)*IM_tbc*E((nx)*(p+1)+1:(nx)*(p+1)+p+1,j);
                SR1_old = SR1_new;
                wh0_old = wh0_new;
                
                u = profile(x,(j-1)*dt);
                U = mxfem1d_FLT( N, FM, u);
                norm_ana = mxfem1d_Norm2( p, N, U);
                eTBC(j)  = mxfem1d_Norm2( p, N, E(:, j)-U);
            end
        end
        
        function r = Gaussian_WP( xx, tt, A0, a, c)
            [x,t] = ndgrid(xx,tt);

            C = (1./(4*a*t*1i+1));
            r = A0.*sqrt(C).*exp(-a*C.*(x-c*t).^2+0.5*1i*c*x-0.25*1i*c^2*t);
        end

        function [y, dy ] = Gaussian_WP1( x, t, A0, a, c)
            % Fixed x
            C  = (1./(4*a*t*1i+1));
            y  = A0.*sqrt(C).*exp(-a*C.*(x-c*t).^2+0.5*1i*c*x-0.25*1i*c^2*t);
            dy = -2*a*(x-c*t).*C.*y+0.5*1i*c*y;
        end
        
        function [err, t] = test_quadW_CQ(c, Nt)
            x_r = 10;
            dx  = 5;
            dt  = 0.5e-3;
            A0  = 1;
            a   = 1/2;
            c   = 8;
            tol = 1e-16;

            t = (0:Nt-1)*dt;

            [y_r, dy_r] = LSE.Gaussian_WP1(x_r, t, A0, a, c);
            [y, dy]     = LSE.Gaussian_WP1(x_r+dx, t, A0, a, c);

            Fdy1 = zeros(1,Nt);
            Fdy2 = zeros(1,Nt);
            Fdy3 = zeros(1,Nt);

            wh1  = LSE.quadW_CQ(Nt, dx, dt, tol);
            wh2  = LSE.quadW_CQ_BDF1(Nt, dx, dt, tol);
            wh3  = LSE.quadW_CQ_BDF2(Nt, dx, dt, tol);

            for i=1:Nt
                wh1_tmp = fliplr(wh1((1:i)));
                wh2_tmp = fliplr(wh2((1:i)));
                wh3_tmp = fliplr(wh3((1:i)));

                y_tmp  = y_r(1:i);
                
                Fdy1(i) = sum(wh1_tmp.*y_tmp);
                Fdy2(i) = sum(wh2_tmp.*y_tmp);
                Fdy3(i) = sum(wh3_tmp.*y_tmp);
            end
            err = zeros(3, Nt);
            err(1,:) = abs(Fdy1 + dy);
            err(2,:) = abs(Fdy2 + dy);
            err(3,:) = abs(Fdy3 + dy);
        end

        function [nx, xi] = fem_index(N, x_l, x_r, x )
            h  = (x_r-x_l)/N;
            d_x = x-x_l;
            nx = floor((x-x_l)/h);
            J  = h/2;
            xm = x_l+nx*h+h/2;
            xi = (x-xm)/J;
        end

        function u = Lobatto(v)
            lv = length(v);
            u(1) = 0.5*(v(1)-v(2));
            u(2) = 0.5*(v(1)+v(2));
    
            for i=3:lv
                u(i) = (1/sqrt(2*(2*(i-1)-1)))*(v(i)-v(i-2));
            end
        end

        function u = diff_legendre(xi,v)
            lv = length(v);
            u(1) = 0;
            for i=2:lv
                u(i) = (i-1)*(v(i-1)-xi*v(i))/(1-xi^2);
            end
        end

        function y = Gaussian_WP3( xx, tt, A0, a, c, g_0, mup)
            [x,t] = ndgrid(xx,tt);

            theta  = mup*t;
            beta   = g_0*sin(theta)/mup;
            nu     = 2.0*g_0*(cos(theta)-1.0)/(mup*mup);
            ibeta2 = g_0*g_0*(theta/2.0 - sin(2*theta)/4.0)/(mup*mup*mup);

            Xi = beta.*x -ibeta2;

            x1 = x+nu-c*t;
            C  = (1./(4*a*t*1i+1));
            y  = A0.*sqrt(C).*exp(-a*C.*(x1).^2+0.5*1i*c*x1+0.25*1i*c^2*t+1i*Xi);
        end
        
        function [y, dy] = Gaussian_WP4( x, t, A0, a, c, g_0, mup)

            theta  = mup*t;
            beta   = g_0*sin(theta)/mup;
            nu     = 2.0*g_0*(cos(theta)-1.0)/(mup*mup);
            ibeta2 = g_0*g_0*(theta/2.0 - sin(2*theta)/4.0)/(mup*mup*mup);

            Xi = beta.*x - ibeta2;
            x1 = x+nu-c*t;

            C  = (1./(4*a*t*1i+1));
            y  = A0.*sqrt(C).*exp(-a*C.*(x1).^2+0.5*1i*c*x1+0.25*1i*c^2*t+1i*Xi);
            dy = -2*a*(x1).*C.*y+0.5*1i*c*y+1i*beta.*y;
            
        end

        function y = V(t, g_0, mup)
            theta = mup*t;
            y     = g_0*cos(theta);
        end

        function y = beta(t, g_0, mup)
            theta = mup*t;
            y     = g_0*sin(theta)/mup;
        end

        function y = nu(t, g_0, mup)
            theta = mup*t;
            y     = 2.0*g_0*(cos(theta)-1.0)/(mup*mup);
        end

        function y = B(t, g_0, mup)
            theta = mup*t;
            y     = g_0*g_0*(theta/2.0 - sin(2*theta)/4.0)/(mup*mup*mup);
        end

        function [err, t] = test_quadW_CQ2(c, Nt)
            x_r = 10;
            dx  = 5;
            x   = x_r + dx;
            dt  = 0.5e-3;
            A0  = 1;
            a   = 1/2;
            % c   = 8;
            tol = 1e-16;

            t = (0:Nt-1)*dt;

            Fdy1  = zeros(1,Nt);
            Fdy2  = zeros(1,Nt);
            Fdy3  = zeros(1,Nt);

            g_0  = 1;
            mup  = 2*pi;
            beta_tmp = LSE.beta(t, g_0, mup);
            nu_tmp   = LSE.nu(t, g_0, mup);
            B_tmp    = LSE.B(t, g_0, mup);
            
            [y_r, dy_r] = LSE.Gaussian_WP4(x_r-nu_tmp, t, A0, a, c, g_0, mup);
            [y, dy]     = LSE.Gaussian_WP4(x, t, A0, a, c, g_0, mup);

            Xi_r = beta_tmp.*(x_r - nu_tmp) - B_tmp;
            Xi  = beta_tmp.*x - B_tmp;
            y2  = y_r.*exp(-1i*Xi_r);
            y3  = (-1i*beta_tmp.*y + dy).*exp(-1i*Xi);

            for i=1:Nt
                wh1  = LSE.quadW_CQ(i, dx+nu_tmp(i), dt, tol);
                wh2  = LSE.quadW_CQ_BDF1(i, dx+nu_tmp(i), dt, tol);
                wh3  = LSE.quadW_CQ_BDF2(i, dx+nu_tmp(i), dt, tol);
                
                wh1_tmp = fliplr(wh1((1:i)));
                wh2_tmp = fliplr(wh2((1:i)));
                wh3_tmp = fliplr(wh3((1:i)));

                y_tmp1  = y2(1:i);

                Fdy1(i)  = sum(wh1_tmp.*y_tmp1);
                Fdy2(i)  = sum(wh2_tmp.*y_tmp1);
                Fdy3(i)  = sum(wh3_tmp.*y_tmp1);
            end
            err = zeros(3, Nt);
            err(1, :) = abs(Fdy1 + y3);
            err(2, :) = abs(Fdy2 + y3);
            err(3, :) = abs(Fdy3 + y3);
        end
    end
end
