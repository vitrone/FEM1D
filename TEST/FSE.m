classdef FSE
    
    properties(Constant)
        CQ_M0 = 20;
    
    end


    methods (Static)
        
        function wh = quadW_CQ(N, dx, dt, tol )
        % Quadrature weights for NtD-map
        % Trapezoidal rule

            if(N<FSE.CQ_M0)
                N = FSE.CQ_M0;
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

        function [E, eTBC] = TBC_CQ_evol(Data)
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

            whr = FSE.quadW_CQ(Nt, dx_r, dt, tol);
            whl = FSE.quadW_CQ(Nt, dx_l, dt, tol);

            eTBC = zeros(1,Nt);
            norm_ana0 = mxfem1d_Norm2( p, N, U);
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
                [nx_l, xi_l] = FSE.fem_index(N, x_l, x_r, x_l+dx_l);
                IMb_l   = mxLGLdataIM( p, xi_l);
                IMvb_l  = FSE.Lobatto(IMb_l);
                GM(1, nx_l+1) = GM(1, nx_l+1) + whl(1)*IMvb_l(1);
                GM(1, nx_l+2) = GM(1, nx_l+2) + whl(1)*IMvb_l(2);
                shift = N+1 + nx_l*(p-2);
                for k=3:p+1
                    GM(1, shift+k) = GM(1, shift+k) + ...
                                       whl(1)*IMvb_l(k);
                end

                [nx_r, xi_r] = FSE.fem_index(N, x_l, x_r, x_r-dx_r);
                IMb_r  = mxLGLdataIM( p, xi_r);
                IMvb_r = FSE.Lobatto(IMb_r);
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

                u = profile(x,(j-1)*dt);
                U = mxfem1d_FLT( N, FM, u);
                % norm_ana = mxfem1d_Norm2( p, N, U);
                eTBC(j)  = mxfem1d_Norm2( p, N, E(:,j)-U)/norm_ana0;
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

            E = U;

            yL = [0];
            yR = [0];

            SL = 0;
            SR = 0;

            whr = FSE.quadW_CQ(Nt, dx_r, dt, tol);
            whl = FSE.quadW_CQ(Nt, dx_l, dt, tol);
            
            eTBC = zeros(1,Nt);

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
                [nx_l, xi_l] = FSE.fem_index(N, x_l, x_r, x_l+dx_l);
                IMb_l   = mxLGLdataIM( p, xi_l);
                IMvb_l  = FSE.Lobatto(IMb_l);
                GM(1, nx_l+1) = GM(1, nx_l+1) + whl(1)*IMvb_l(1);
                GM(1, nx_l+2) = GM(1, nx_l+2) + whl(1)*IMvb_l(2);
                shift = N+1 + nx_l*(p-2);
                for k=3:p+1
                    GM(1, shift+k) = GM(1, shift+k) + ...
                                       whl(1)*IMvb_l(k);
                end

                [nx_r, xi_r] = FSE.fem_index(N, x_l, x_r, x_r-dx_r);
                IMb_r  = mxLGLdataIM( p, xi_r);
                IMvb_r = FSE.Lobatto(IMb_r);
                GM(N+1, nx_r+1) = GM(N+1, nx_r+1) + whr(1)*IMvb_r(1);
                GM(N+1, nx_r+2) = GM(N+1, nx_r+2) + whr(1)*IMvb_r(2);
                shift = N+1 + nx_r*(p-2);
                for k=3:p+1
                    GM(N+1, shift+k) = GM(N+1, shift+k) + ...
                                       whr(1)*IMvb_r(k);
                end
            end
            
            norm_ana0 = mxfem1d_Norm2( p, N, U);
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
                % norm_ana = mxfem1d_Norm2( p, N, U);
                eTBC(j)  = mxfem1d_Norm2( p, N, E-U)/norm_ana0;
            end
        end

        function r = Gaussian_WP( xx, tt, A0, a, c)
            [x,t] = ndgrid(xx,tt);

            C = (1./(4*a*t*1i+1));
            r = A0.*sqrt(C).*exp(-a*C.*(x-c*t).^2+0.5*1i*c*x-0.25*1i*c^2*t);
        end

        % Hermite-Gauss wavepacket Solutions
        function r = HG2_WP( xx, tt, A0, a, c)
            [x,t] = ndgrid(xx,tt);

            C = (1./(4*a*t*1i+1));
            x1 = x-c*t;
            hpoly = 4*x1.^2 -2;
            theta = atan(4*a*t);
            E = exp(-a*C.*(x-c*t).^2+0.5*1i*c*x-0.25*1i*c^2*t)
            r = (A0/(4*a^(1/2)*pi^(1/4)))*exp(-1i*2*theta).*hpoly.*sqrt(C).*E;
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

    end
end
