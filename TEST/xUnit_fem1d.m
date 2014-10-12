classdef xUnit_fem1d < matlab.unittest.TestCase
% test = test_fem1d
% test.run
% run(test, 'test_GP_ABC1a_CQ(testCase')

    methods(Test)
        function test_interpolation(testCase)
            x_l = -5;
            x_r =  5;
            N   =  10;
            p   =  4;
            poly_fun = @(x) (1+1i)*x.^p;
            tol = 1e-9;
            [xir, quadW, IM, FM] = mxLGLdataLT2( p, tol);
            x = ref2mesh (xir, N, [x_l x_r]);
            
            u = poly_fun(x);      % a column vector
            U = mxfem1d_FLT( N, FM, u);

            xi = (-1:0.001:1)';
            P  = length(xi)-1;
            IM1 = mxLGLdataIM( p, xi);
            u1  = mxfem1d_ILT( N, IM1, U);
            x1  = ref2mesh (xi, N, [x_l x_r]);
            u_actual = poly_fun(x1);
            e = norm(u1-u_actual)/norm(u_actual);
            % disp(e)
            testCase.verifyLessThan( e, tol );
        end

        function test_interpolation2(testCase)
            x_l = -5;
            x_r =  5;
            N   =  10;
            p   =  4;
            poly_fun = @(x) (1+1i)*[x.^(p-2), x.^(p-1), x.^(p)];
            tol = 1e-9;
            [xir, quadW, IM, FM] = mxLGLdataLT2( p, tol);
            x = ref2mesh (xir, N, [x_l x_r]);
            
            u = poly_fun(x);      % a column vector
            U = mxfem1d_FLT( N, FM, u);

            xi = (-1:0.001:1)';
            P  = length(xi)-1;
            IM1 = mxLGLdataIM( p, xi);
            u1  = mxfem1d_ILT( N, IM1, U);
            x1  = ref2mesh (xi, N, [x_l x_r]);
            u_actual = poly_fun(x1);
            e(1) = norm(u1(:,1)-u_actual(:,1))/norm(u_actual(:,1));
            e(2) = norm(u1(:,2)-u_actual(:,2))/norm(u_actual(:,2));
            e(3) = norm(u1(:,3)-u_actual(:,3))/norm(u_actual(:,3));
            % e = norm(u1-u_actual)./norm(u_actual);
            % disp(e)
            testCase.verifyLessThan( e, tol );
        end

        function test_norm(testCase)
            x_l = -5;
            x_r =  5;
            N   =  20;
            J   =  (x_r-x_l)/(2*N);
            p   =  20;
            fun = @(x) exp(-x.^2/2);
            expected_norm = (pi)^(1/4);
            tol = 1e-9;
            [xir, quadW, IM, FM] = mxLGLdataLT2( p, tol);
            x = ref2mesh (xir, N, [x_l x_r]);
            
            u = fun(x);      % a column vector
            U = mxfem1d_FLT( N, FM, u);
            lp_norm = sqrt(J)*mxfem1d_Norm2(p, N, U);

            e = abs(expected_norm-lp_norm)/abs(expected_norm);
            % disp(e)
            testCase.verifyLessThan( e, tol );
        end
        
        function test2_norm(testCase)
            x_l = -5;
            x_r =  5;
            N   =  20;
            J   =  (x_r-x_l)/(2*N);
            p   =  8;
            fun = @(x) exp(-x.^2/2);
            expected_norm = (pi)^(1/4);
            tol = 1e-9;
            [xir, quadW, IM, FM] = mxLGLdataLT2( p, tol);
            x = ref2mesh (xir, N, [x_l x_r]);
            
            u = fun(x);      % a column vector
            U = mxfem1d_FLT( N, FM, u);
            lp_norm = sqrt(J)*mxfem1d_Norm2(p, N, U);

            e = abs(expected_norm-lp_norm)/abs(expected_norm);
            % disp(e)
            testCase.verifyLessThan( e, tol );
        end

        function test_fem1d_L2F(testCase)
            x_l = -5;
            x_r =  5;
            N   =  40;

            fun = @(x) exp(-(1+1i)*x.^2/2);
            tol = 1e-9;

            for p=2:12
                [xi, quadW] = mxLGLdataLT1( p, tol);

                FM = mxLGLdataFM(p, xi);
                IM = mxLGLdataIM(p, xi);
                x  = ref2mesh(xi, N, [x_l x_r]);

                u = fun(x);      
                U1 = mxfem1d_FLT( N, FM, u);
                vb = mxfem1d_L2F(p, U1);
                U2 = mxfem1d_F2L(p, vb);
                
                e = norm(U2-U1)/norm(U1);
                testCase.verifyLessThan( e, tol );
                disp(e)
            end

        end
        function test_fem1d_sparse_GMM(testCase)
            x_l = -5;
            x_r =  5;
            N   =  600;
            p   =  4;
            nr_LGL = 2*p+1;
            P = nr_LGL-1;

            fun = @(x) exp(-(1+1i)*x.^2/2);
            pot = @(x) 1+1i*(x);
            tol = 1e-9;

            [xi, quadW] = mxLGLdataLT1( P, tol);

            FM = mxLGLdataFM(p, xi);
            IM = mxLGLdataIM(p, xi);
            Q  = mxfem1d_quadM(quadW, IM);
            x  = ref2mesh (xi, N, [x_l x_r]);

            % u -> a column vector
            u = fun(x);      
            U = mxfem1d_FLT( N, FM, u);
            vb = mxfem1d_L2F(p, U);

            phi = pot(x);
            M = mxfem1d_sparse_GMM(p, Q, phi);

            u1 = u.*phi;
            U1 = mxfem1d_FLT( N, FM, u1);
            Pvb = mxfem1d_PrjL2F(p, U1);

            sM = M + transpose(M) - diag(diag(M));
            Pvb1 = sM*vb;
            % disp(sM)
            e = norm(Pvb-Pvb1)/norm(Pvb);
            testCase.verifyLessThan( e, tol );
            disp(e)
        end

        function test_GP_ABC1a_CQ(testCase)
            I = sqrt(-1);
            tol = 1e-9; % tolerance for LGL points
            N   = 400;  % Number of finite elements
            Nth = 50;
            p   = 4;    % highest degree of polynomials
            [xr quadW IM FM] = mxLGLdataLT2( p, tol);
            % xr: LGL points
            % FM: matrix for forward transform
            % IM: matrix for inverse transform

            P = 2*p; % For the purpose of numerical integration
            % Degree of the Legendre polynomials involved --> P
            [xri quadWi] = mxLGLdataLT1( P, tol);
            FMi = mxLGLdataFM( p, xri); 
            IMi = mxLGLdataIM( p, xri);

            dt  =  0.5e-3;
            r   =  2/dt;
            Nt  =  2000+1;
            x_r =  10;
            x_l = -10;
            x   =  ref2mesh(xr,N,[x_l x_r]);
            
            a0  =  1;
            eta =  1;
            c   =  15;
            chi =  2;
            KL  =  c/2;
            KR  =  KL;
            AL  =  a0*exp(I*KL*x_l);
            AR  =  a0*exp(I*KR*x_r);

            mu  = 2*pi;
            g_0 = 0;
            g_1 = 1;
            pot = sprintf('@(t)%d+%d*cos(%d*t)', g_0,g_1,mu);
            profile = sprintf( '@(x,t)GP_EQUATIONS.gpe_bs(x,t,%d,%d,%d,%d,%d,%d)',...
                               a0, eta, c, g_0, g_1,mu);
            hXiL    = sprintf( '@(t) GP_EQUATIONS.Xi(%d,t,%d,%d,%d,%d,%d,%d)',...
                               x_l,a0,KL,g_0,g_1,mu,chi);
            hXiR    = sprintf( '@(t) GP_EQUATIONS.Xi(%d,t,%d,%d,%d,%d,%d,%d)',...
                               x_r,a0,KR,g_0,g_1,mu,chi);
            hdXi    = sprintf('@(t) GP_EQUATIONS.dXi(t,%d,%d,%d)',g_0,g_1,mu);
            Data = struct( 'fN',N,...
                           'fp',p,...
                           'fP',P,...
                           'fFM',FM,...
                           'fIM',IM,...
                           'fFMi',FMi,...
                           'fIMi',IMi,...
                           'fquadWi',quadWi,...
                           'fxri',xri,...
                           'fdt',dt,...
                           'fNt',Nt,...
                           'fNth',Nth,...
                           'fx_r',x_r,...
                           'fx_l',x_l,...
                           'fx',x,...
                           'fAL',AL,...
                           'fAR',AR,...
                           'fKL',KL,...
                           'fKR',KR,...
                           'fchi',chi,...
                           'hpot',pot,...
                           'hXiL',hXiL,...
                           'hXiR',hXiR,...
                           'hdXi',hdXi,...
                           'hfunc',profile,...
                           'ftol',tol);

            
            E = GP_EQUATIONS.ABC1a_CQ1(Data);
            e = max(E);
            disp(e)
            testCase.verifyLessThan( e, 3e-3 );
        end

        function test_GP_ABC1a_Pade(testCase)
            I = sqrt(-1);
            tol = 1e-9; % tolerance for LGL points
            N   = 400;  % Number of finite elements
            Nth = 50;
            p   = 4;    % highest degree of polynomials
            [xr quadW IM FM] = mxLGLdataLT2( p, tol);
            % xr: LGL points
            % FM: matrix for forward transform
            % IM: matrix for inverse transform

            P = 2*p; % For the purpose of numerical integration
            % Degree of the Legendre polynomials involved --> P
            [xri quadWi] = mxLGLdataLT1( P, tol);
            FMi = mxLGLdataFM( p, xri); 
            IMi = mxLGLdataIM( p, xri);

            dt  =  0.5e-3;
            r   =  2/dt;
            Nt  =  2000+1;
            x_r =  10;
            x_l = -10;
            x   =  ref2mesh(xr,N,[x_l x_r]);
            Mp  =  50; 
            
            a0  =  1;
            eta =  1;
            c   =  15;
            chi =  2;
            KL  =  c/2;
            KR  =  KL;
            AL  =  a0*exp(I*KL*x_l);
            AR  =  a0*exp(I*KR*x_r);

            mu  = 2*pi;
            g_0 = 0;
            g_1 = 1;
            pot = sprintf('@(t)%d+%d*cos(%d*t)', g_0,g_1,mu);
            profile = sprintf( '@(x,t)GP_EQUATIONS.gpe_bs(x,t,%d,%d,%d,%d,%d,%d)',...
                               a0, eta, c, g_0, g_1,mu);
            hXiL    = sprintf( '@(t) GP_EQUATIONS.Xi(%d,t,%d,%d,%d,%d,%d,%d)',...
                               x_l,a0,KL,g_0,g_1,mu,chi);
            hXiR    = sprintf( '@(t) GP_EQUATIONS.Xi(%d,t,%d,%d,%d,%d,%d,%d)',...
                               x_r,a0,KR,g_0,g_1,mu,chi);
            hdXi    = sprintf('@(t) GP_EQUATIONS.dXi(t,%d,%d,%d)',g_0,g_1,mu);
            Data = struct( 'fN',N,...
                           'fp',p,...
                           'fP',P,...
                           'fFM',FM,...
                           'fIM',IM,...
                           'fFMi',FMi,...
                           'fIMi',IMi,...
                           'fquadWi',quadWi,...
                           'fxri',xri,...
                           'fdt',dt,...
                           'fNt',Nt,...
                           'fNth',Nth,...
                           'fMp',Mp,...
                           'fx_r',x_r,...
                           'fx_l',x_l,...
                           'fx',x,...
                           'fAL',AL,...
                           'fAR',AR,...
                           'fKL',KL,...
                           'fKR',KR,...
                           'fchi',chi,...
                           'hpot',pot,...
                           'hXiL',hXiL,...
                           'hXiR',hXiR,...
                           'hdXi',hdXi,...
                           'hfunc',profile,...
                           'ftol',tol);

            
            E = GP_EQUATIONS.ABC1a_Pade1(Data);
            e = max(E);
            disp(e)
            testCase.verifyLessThan( e, 3e-3 );
        end
    end
end

