classdef LSE
   
% Solving Linear Schroedinger equation using FEM
%   iu_t + u_xx + phi(x)u + V(t)xu = 0.
%
% Boundary Conditions:
%   Natural boundary condition
%   Transparent Boundary conditons
% This classes uses another class of objects LSE_WP which defined 
% wave-packet solutions (WP) for the LSE considered above.
% See LSE_WP
%
% Properties
%   |--tol             : Tolerance
%   |--p               : Highest degree of the legendre polynomial
%   |--P               : Polynomial degree for quadrature 
%   |--N               : Number of finite-elements
%   |--Nt              : Number of time-steps
%   |--t               : 1-by-Nt vector storing time-grid points
%   |--x               : Space grid
%   |--dt              : Size of the time-step
%   |--domain          : Computational dimain [x_l, x_r]
%   |--xpotential      : Static potential
%   |--lin_tpotential  : Struct, time-dependent potential 
%   |--E               : Field in Legendre basis
%   |--eTBC            : Error
%   |--soln_name       
%   |--soln_params
%   |--WP_data
%
% EXAMPLE:
%    tol = 1e-9; % tolerance for LGL points
%    N   = 400;  % Number of finite elements
%    p   = 4;    % highest degree of polynomials
%    
%    dt  =  0.5e-3;
%    Nt  =  4000;
%    obj = LSE(p, N, dt, Nt, tol);
%    
%    obj.domain = [-10, 10];
%    t0    = 1;
%    sigma = 1/8;
%    tpotential = struct( 'g_0'            , 1,...
%                         'mup'            , 2*pi,...
%                         't_0'            , t0,...
%                         'sigma'          , sigma,... 
%                         'potential_name' , 'GAUSSIAN',...
%                         'T_supp'         , t0 + 8*sigma);
%
%    obj.lin_tpotential = tpotential;
%    obj = obj.init_solver();
%    soln_params = struct( 'A0', 1,...
%                          'c' , 8,...
%                          'a' , 1);
%    
%    options = struct('BC_type'     , 'TBC_CQ',...
%                     'mode'        , 'ERROR',...
%                     'soln_name'   , 'Gaussian_WP',...
%                     'soln_params' , soln_params);
%    obj = obj.solve(options);
%    logerr_plot(obj.t, obj.eTBC)

    properties(SetAccess=private)
        tol = 1e-9;
        p
        P
        N
        Nt
        t
        x
        dt
    end
    properties
        domain         % computational domain: [x_l, x_r]
        xpotential     % Static potential
        lin_tpotential % time-dependent linear potential
        E              % Field in Legendre basis
        E_ana          % Field in Legendre basis, analytical
        eTBC           % L2 norm of error
        soln_name      % name of the analytical solution
        soln_params    % parameters 
        WP_data
    end
    properties(SetAccess = private, GetAccess = private)
        u0
        dx
        J   % Jacobian of transformation
        r
        sw % defines the strip width  
        % xr   : LGL points
        % quadW: Quadrature weights for LGL points
        % FM   : matrix for forward transform
        % IM   : matrix for inverse transform
        LGLdata
        % xri   : LGL points
        % quadWi: Quadrature weights for LGL points
        % FMi   : matrix for forward transform
        % IMi   : matrix for inverse transform
        LGLdataQ
        x1
        GM
        GPMM
        phi
        V

        nr_fem
        nr_lgl

        tbeta
        tnu
        tB
    end

    properties(Constant, SetAccess = private, GetAccess = private)
        CQ_M0 = 20;
        stol = 1e-4;
    
    end

    methods
        function obj = LSE(p, N, dt, Nt, tol)
            obj.p = p;
            obj.N = N;
            obj.tol = tol;
            [obj.LGLdata.xr, obj.LGLdata.quadW, obj.LGLdata.IM, obj.LGLdata.FM] ...
                                                            = mxLGLdataLT2( p, tol);
            % xr: LGL points
            % FM: matrix for forward transform
            % IM: matrix for inverse transform

            obj.P = 2*p; % For the purpose of numerical integration
            % Degree of the Legendre polynomials involved --> P
            [obj.LGLdataQ.xri, obj.LGLdataQ.quadWi] = mxLGLdataLT1( obj.P, tol);
            obj.LGLdataQ.FMi = mxLGLdataFM( p, obj.LGLdataQ.xri);
            obj.LGLdataQ.IMi = mxLGLdataIM( p, obj.LGLdataQ.xri);
            
            obj.dt = dt;
            obj.Nt = Nt;
            obj.t  = ((0:Nt-1)*dt);
        end

        function obj = set.domain(obj, domain)
            obj.domain = domain;
            obj.x   = ref2mesh(obj.LGLdata.xr  , obj.N, domain);
            obj.x1  = ref2mesh(obj.LGLdataQ.xri, obj.N, domain);
            obj.phi = [];
            obj.dx  = (obj.domain(2)-obj.domain(1))/obj.N;
            obj.J   = obj.dx/2;
            obj.r   = 2/obj.dt;
        end

        function obj = set.lin_tpotential(obj, tpotential)
        % Define the laser field
        % USER_DEFINED - require function handle 
        % CONSTANT     - constant potential
        % Gaussian     - Gaussian pulse, finite support
        % CW           - Continuous Wave

            obj.lin_tpotential = tpotential;
            switch(tpotential.potential_name)
                case{'USER_DEFINED'}
                    obj.V = tpotential.func(obj.t);
                case{'CONSTANT'}
                    obj.V = ones(1, obj.Nt);
                case{'CW'}
                    obj.V = LSE.tpotential_CW(obj.t, tpotential);
                case{'GAUSSIAN'}
                    obj.V = LSE.tpotential_Gaussian(obj.t, tpotential);
                otherwise
                    error('Unknown potential')
            end

            %obj.WP_data = LSE_WP(obj.dt, obj.Nt, params);
            if isinf(tpotential.T_supp)
                [obj.tnu, obj.tbeta, obj.tB] = LSE.auxilary_func( obj.dt,...
                                                                  obj.V);
            else
                [obj.tnu, obj.tbeta, obj.tB] = LSE.auxilary_func( obj.dt,...
                                                                  obj.V,...
                                                                  tpotential.T_supp);
            end

            K2 =  max(0, max(obj.tnu)) + 1.03;
            K1 = -min(0, min(obj.tnu)) + 1.03;
            obj.sw = [K2, K1];
        end

        function obj = set.xpotential(obj, xpotential)
            obj.xpotential = xpotential;
            obj.phi = xpotential(obj.x1);
        end

        function obj = init_solver(obj)
        % Initialize the solver 

            Q = mxfem1d_quadM( obj.LGLdataQ.quadWi, obj.LGLdataQ.IMi);
            pot = ones(size(obj.x1));
            % mass matrix
            LGMM = mxfem1d_sparse_GMM( obj.p, Q, pot);
            UGMM = transpose(LGMM) - diag(diag(LGMM));
            
            %% stiffness matrix
            [row1, col1, gsm, nnz1_] = fem1d_crgsm( obj.p, obj.N);
            GSM = sparse(row1, col1, gsm, obj.N*obj.p+1, obj.N*obj.p+1, nnz1_);

            LGPMM = mxfem1d_sparse_GMM( obj.p, Q, obj.x1);
            UGPMM = transpose(LGPMM) - diag(diag(LGPMM));
            obj.GPMM = -obj.J*(UGPMM+LGPMM);  

            if ~isempty(obj.phi)
                LGPMM = mxfem1d_sparse_GMM( obj.p, Q, obj.phi);
                UGPMM = transpose(LGPMM) - diag(diag(LGPMM));
                GPMMs  = -obj.J*(UGPMM+LGPMM); 
                obj.GM = GSM/obj.J + (-1i*obj.r*obj.J)*(UGMM+LGMM) + GPMMs;
            else
                obj.GM = GSM/obj.J + (-1i*obj.r*obj.J)*(UGMM+LGMM);
            end
            %%
            obj.nr_fem = obj.p*obj.N + 1;
            obj.nr_lgl = (obj.p + 1)*obj.N;
        end

        function obj = solve(obj, options)
            % xx
            % 11: EVOLVE, NATURAL
            % 21: ERROR , NATURAL
            % 12: EVOLVE, TBC_CQ
            % 22: ERROR , TBC_CQ
            action_code = 00; 

            switch(options.BC_type)
                case{'NATURAL'}
                    action_code = action_code + 1;
                case{'TBC_CQ'}
                    action_code = action_code + 2;
                otherwise
                    error('Boundary condition type is unspecified!')
            end
            switch(options.mode)
                case{'EVOLVE'}
                    action_code = action_code + 10;
                    obj.u0 = options.u0;
                case{'ERROR'}
                    action_code = action_code + 20;
                    if(isinf(obj.lin_tpotential.T_supp))
                        obj.WP_data = LSE_WP( obj.dt,...
                                              obj.Nt,...
                                              obj.lin_tpotential,...
                                              obj.V);
                    else
                        obj.WP_data = LSE_WP( obj.dt,...
                                              obj.Nt,...
                                              obj.lin_tpotential);
                    end
                    obj.soln_name   = options.soln_name;
                    obj.soln_params = options.soln_params;
                case{'EVOLVE_ERROR'}
                    action_code = action_code + 30;
                    if(isinf(obj.lin_tpotential.T_supp))
                        obj.WP_data = LSE_WP( obj.dt,...
                                              obj.Nt,...
                                              obj.lin_tpotential,...
                                              obj.V);
                    else
                        obj.WP_data = LSE_WP( obj.dt,...
                                              obj.Nt,...
                                              obj.lin_tpotential);
                    end
                    obj.soln_name   = options.soln_name;
                    obj.soln_params = options.soln_params;
                otherwise
                    error('Solver mode is unspecified!')
            end
            switch(action_code)
                case{11}
                    obj = obj.Natural_evolve();
                case{21}
                    obj = obj.Natural_error();
                case{12}
                    obj = obj.TBC_CQ_evolve();
                case{22}
                    obj = obj.TBC_CQ_error();
                case{32}
                    obj = obj.TBC_CQ_evolve_error();
            end
        end

        function plot_error(obj)
            if ~isempty(obj.eTBC)
                logerr_plot(obj.t, obj.eTBC)
            end
        end

    end
    
    methods(Access = private)

        function obj = Natural_evolve(obj)
            % Solving Linear Schroedinger equation using FEM
            % iu_t+u_xx+V(t)xu=0
            % Natural boundary condition
            % Argument structure is 
            %
            % out ---> N x Nt matrix containing Legendre transform of the solution on 
            % the FEM grid.
            %

            U  = mxfem1d_FLT( obj.N, obj.LGLdata.FM, obj.u0);
            X  = zeros(obj.nr_fem, 1);
            Y  = zeros(obj.nr_lgl, 1);

            obj.E = zeros(obj.nr_lgl, obj.Nt);
            obj.E(:, 1) = U;

            for j=2:obj.Nt,

                tmpV = 0.5*(obj.V(j-1) + obj.V(j));
                GM1 = obj.GPMM*tmpV + obj.GM;

                Pvb = mxfem1d_PrjL2F( obj.p, obj.E(:, j-1));
                F   = -1i*obj.r*obj.J*Pvb;

                X = GM1\F;
                Y = mxfem1d_F2L( obj.p, X);
                obj.E(:, j)  = 2*Y-obj.E(:, j-1);
            end
        end

        function obj = Natural_error(obj)
            % Solving Linear Schroedinger equation using FEM
            % iu_t+u_xx+V(t)xu=0
            % Natural boundary condition

            u = obj.WP_data.(obj.soln_name)(obj.x, 1, obj.soln_params);
            U = mxfem1d_FLT( obj.N, obj.LGLdata.FM, u);
            X = zeros(obj.nr_fem, 1);
            Y = zeros(obj.nr_lgl, 1);

            E = zeros(obj.nr_lgl, 1);
            E(:,1) = U;
            obj.eTBC = zeros(1, obj.Nt);
            norm0 = mxfem1d_Norm2( obj.p, obj.N, U);
            for j=2:obj.Nt,

                tmpV = 0.5*(obj.V(j-1) + obj.V(j));
                GM1 = obj.GPMM*tmpV + obj.GM;

                Pvb = mxfem1d_PrjL2F( obj.p, E);
                F   = -1i*obj.r*obj.J*Pvb;

                X = GM1\F;
                Y = mxfem1d_F2L( obj.p, X);
                E  = 2*Y-E;

                u = obj.WP_data.(obj.soln_name)(obj.x, j, obj.soln_params);
                U = mxfem1d_FLT( obj.N, obj.LGLdata.FM, u);
                obj.eTBC(j) = mxfem1d_Norm2( obj.p, obj.N, E-U)/norm0;
            end
        end

        function obj = TBC_CQ_evolve(obj)
            dx_l    = obj.sw(1);
            dx_r    = obj.sw(2);

            % verify support of the initial data
            supp = LSE.find_support(obj.x, obj.u0, obj.stol);
            if ((obj.domain(1)>supp(1)-dx_l) || (obj.domain(2)<supp(2)+dx_r))
                fprintf( 'Support of the initial data: [%0.5f, %0.5f]\n',...
                          supp(1), supp(2));
                error('Computational domain is inadmissible!')
            end
            if~isempty(obj.phi)
                supp = LSE.find_support(obj.x1, obj.phi, obj.stol);
                if ((obj.domain(1)>supp(1)-dx_l) || (obj.domain(2)<supp(2)+dx_r))
                    fprintf( 'Support of the static potential: [%0.5f, %0.5f]\n',...
                              supp(1), supp(2));
                    error('Computational domain is inadmissible!')
                end
            end
            U  = mxfem1d_FLT( obj.N, obj.LGLdata.FM, obj.u0);
            X  = zeros(obj.nr_fem, 1);
            Y  = zeros(obj.nr_lgl, 1);

            obj.E = zeros(obj.nr_lgl, obj.Nt);
            obj.E(:,1) = U;

            xl  = obj.domain(1) + dx_l - obj.tnu;
            xr  = obj.domain(2) - dx_r - obj.tnu;

            yL = [0];
            yR = [0];

            SL = 0;
            SR = 0;

            
            whl0_old = 0;
            whr0_old = 0;

            SR_old  = 0;
            SL_old  = 0;
            for j=2:obj.Nt,

                nu_tmp   = 0.5*(obj.tnu(j)   + obj.tnu(j-1));
                beta_tmp = 0.5*(obj.tbeta(j) + obj.tbeta(j-1));
                B_tmp    = 0.5*(obj.tB(j)    + obj.tB(j-1));

                whl = LSE.quadW_CQ(j, dx_l - obj.tnu(j), obj.dt, obj.tol);
                whr = LSE.quadW_CQ(j, dx_r + obj.tnu(j), obj.dt, obj.tol);

                tmpV = 0.5*(obj.V(j-1) + obj.V(j));
                GM1 = obj.GPMM*tmpV + obj.GM;

                xl_tmp  = 0.5*(xl(j) + xl(j-1));
                xr_tmp  = 0.5*(xr(j) + xr(j-1));

                [nx_l, xi_l] = LSE.fem_index(obj.N, obj.domain, xl_tmp);
                [nx_r, xi_r] = LSE.fem_index(obj.N, obj.domain, xr_tmp);

                IMb_l = mxLGLdataIM( obj.p, xi_l);
                IMb_r = mxLGLdataIM( obj.p, xi_r);
                 
                whl0_new = whl(1)*exp(1i*obj.tbeta(j)*( -dx_l + obj.tnu(j)));
                whr0_new = whr(1)*exp(1i*obj.tbeta(j)*(  dx_r + obj.tnu(j)));

                IMvb_l = 0.5*(whl0_new+whl0_old)*LSE.Lobatto(IMb_l);
                IMvb_r = 0.5*(whr0_new+whr0_old)*LSE.Lobatto(IMb_r);

                GM1(1, 1)      = GM1(1, 1)  + 1i*beta_tmp;
                GM1(1, nx_l+1) = GM1(1, nx_l+1) + IMvb_l(1);
                GM1(1, nx_l+2) = GM1(1, nx_l+2) + IMvb_l(2);
                
                shift = obj.N+1 + nx_l*(obj.p-2);

                for k=3:obj.p+1
                    GM1(1, shift+k) = GM1( 1, shift+k) + ...
                                           IMvb_l(k);
                end

                GM1(obj.N+1, obj.N+1) = GM1(obj.N+1, obj.N+1)  - 1i*beta_tmp;
                GM1(obj.N+1, nx_r+1)  = GM1(obj.N+1, nx_r+1) + IMvb_r(1);
                GM1(obj.N+1, nx_r+2)  = GM1(obj.N+1, nx_r+2) + IMvb_r(2);
                
                shift = obj.N + 1 + nx_r*(obj.p-2);
                
                for k=3:obj.p+1
                    GM1(obj.N+1, shift+k) = GM1( obj.N+1, shift+k) + ...
                                             IMvb_r(k);
                end

                Pvb = mxfem1d_PrjL2F( obj.p, obj.E(:, j-1));
                F   = -1i*obj.r*obj.J*Pvb;


                whl_tmp = fliplr(whl((2:j)));
                whr_tmp = fliplr(whr((2:j)));

                Xi1_l = obj.tbeta(j)*xl(j)         - obj.tB(j); 
                Xi_l  = obj.tbeta(j)*obj.domain(1) - obj.tB(j); 

                Xi1_r = obj.tbeta(j)*xr(j)         - obj.tB(j); 
                Xi_r  = obj.tbeta(j)*obj.domain(2) - obj.tB(j); 

                SL_new = exp(1i*Xi_l)*sum(whl_tmp.*yL);
                SR_new = exp(1i*Xi_r)*sum(whr_tmp.*yR);


                tmpF = F;
                tmpF(1)       = F(1)       - 0.5*(SL_new + SL_old);
                tmpF(obj.N+1) = F(obj.N+1) - 0.5*(SR_new + SR_old);

                X = GM1\tmpF;
                Y = mxfem1d_F2L( obj.p, X);
                obj.E(:, j) = 2*Y-obj.E(:, j-1);

                indl = (nx_l)*(obj.p+1)+1:(nx_l)*(obj.p+1)+obj.p+1;
                indr = (nx_r)*(obj.p+1)+1:(nx_r)*(obj.p+1)+obj.p+1;
                yL(1, j) = exp(-1i*Xi1_l)*IMb_l*obj.E(indl,j);
                yR(1, j) = exp(-1i*Xi1_r)*IMb_r*obj.E(indr,j);
                
                SR_old = SR_new;
                SL_old = SL_new;
                
                whl0_old = whl0_new;
                whr0_old = whr0_new;
            end
        end


        function obj = TBC_CQ_evolve_error(obj)
            dx_l    = obj.sw(1);
            dx_r    = obj.sw(2);

            % verify support of the initial data
            u = obj.WP_data.(obj.soln_name)(obj.x, 1, obj.soln_params);

            supp = LSE.find_support(obj.x, u, obj.stol);
            
            if ((obj.domain(1)>supp(1)-dx_l) || (obj.domain(2)<supp(2)+dx_r))
                fprintf( 'Support of the initial data: [%0.5f, %0.5f]\n',...
                          supp(1), supp(2));
                error('Computational domain is inadmissible!')
            end
            if~isempty(obj.phi)
                supp = LSE.find_support(obj.x1, obj.phi, obj.stol);
                if ((obj.domain(1)>supp(1)-dx_l) || (obj.domain(2)<supp(2)+dx_r))
                    fprintf( 'Support of the static potential: [%0.5f, %0.5f]\n',...
                              supp(1), supp(2));
                    error('Computational domain is inadmissible!')
                end
            end
            U  = mxfem1d_FLT( obj.N, obj.LGLdata.FM, u);
            X  = zeros(obj.nr_fem, 1);
            Y  = zeros(obj.nr_lgl, 1);

            obj.E     = zeros(obj.nr_lgl, obj.Nt);
            obj.E_ana = zeros(obj.nr_lgl, obj.Nt);
            obj.E(:, 1)     = U;
            obj.E_ana(:, 1) = U;

            obj.eTBC = zeros(1, obj.Nt);
            norm0 = mxfem1d_Norm2( obj.p, obj.N, U);

            xl  = obj.domain(1) + dx_l - obj.tnu;
            xr  = obj.domain(2) - dx_r - obj.tnu;

            yL = [0];
            yR = [0];

            SL = 0;
            SR = 0;

            
            whl0_old = 0;
            whr0_old = 0;

            SR_old  = 0;
            SL_old  = 0;
            for j=2:obj.Nt,

                nu_tmp   = 0.5*(obj.tnu(j)   + obj.tnu(j-1));
                beta_tmp = 0.5*(obj.tbeta(j) + obj.tbeta(j-1));
                B_tmp    = 0.5*(obj.tB(j)    + obj.tB(j-1));

                whl = LSE.quadW_CQ(j, dx_l - obj.tnu(j), obj.dt, obj.tol);
                whr = LSE.quadW_CQ(j, dx_r + obj.tnu(j), obj.dt, obj.tol);

                tmpV = 0.5*(obj.V(j-1) + obj.V(j));
                GM1 = obj.GPMM*tmpV + obj.GM;

                xl_tmp  = 0.5*(xl(j) + xl(j-1));
                xr_tmp  = 0.5*(xr(j) + xr(j-1));

                [nx_l, xi_l] = LSE.fem_index(obj.N, obj.domain, xl_tmp);
                [nx_r, xi_r] = LSE.fem_index(obj.N, obj.domain, xr_tmp);

                IMb_l = mxLGLdataIM( obj.p, xi_l);
                IMb_r = mxLGLdataIM( obj.p, xi_r);
                 
                whl0_new = whl(1)*exp(1i*obj.tbeta(j)*( -dx_l + obj.tnu(j)));
                whr0_new = whr(1)*exp(1i*obj.tbeta(j)*(  dx_r + obj.tnu(j)));

                IMvb_l = 0.5*(whl0_new+whl0_old)*LSE.Lobatto(IMb_l);
                IMvb_r = 0.5*(whr0_new+whr0_old)*LSE.Lobatto(IMb_r);

                GM1(1, 1)      = GM1(1, 1)  + 1i*beta_tmp;
                GM1(1, nx_l+1) = GM1(1, nx_l+1) + IMvb_l(1);
                GM1(1, nx_l+2) = GM1(1, nx_l+2) + IMvb_l(2);
                
                shift = obj.N+1 + nx_l*(obj.p-2);

                for k=3:obj.p+1
                    GM1(1, shift+k) = GM1( 1, shift+k) + ...
                                           IMvb_l(k);
                end

                GM1(obj.N+1, obj.N+1) = GM1(obj.N+1, obj.N+1)  - 1i*beta_tmp;
                GM1(obj.N+1, nx_r+1)  = GM1(obj.N+1, nx_r+1) + IMvb_r(1);
                GM1(obj.N+1, nx_r+2)  = GM1(obj.N+1, nx_r+2) + IMvb_r(2);
                
                shift = obj.N + 1 + nx_r*(obj.p-2);
                
                for k=3:obj.p+1
                    GM1(obj.N+1, shift+k) = GM1( obj.N+1, shift+k) + ...
                                             IMvb_r(k);
                end

                Pvb = mxfem1d_PrjL2F( obj.p, obj.E(:, j-1));
                F   = -1i*obj.r*obj.J*Pvb;


                whl_tmp = fliplr(whl((2:j)));
                whr_tmp = fliplr(whr((2:j)));

                Xi1_l = obj.tbeta(j)*xl(j)         - obj.tB(j); 
                Xi_l  = obj.tbeta(j)*obj.domain(1) - obj.tB(j); 

                Xi1_r = obj.tbeta(j)*xr(j)         - obj.tB(j); 
                Xi_r  = obj.tbeta(j)*obj.domain(2) - obj.tB(j); 

                SL_new = exp(1i*Xi_l)*sum(whl_tmp.*yL);
                SR_new = exp(1i*Xi_r)*sum(whr_tmp.*yR);


                tmpF = F;
                tmpF(1)       = F(1)       - 0.5*(SL_new + SL_old);
                tmpF(obj.N+1) = F(obj.N+1) - 0.5*(SR_new + SR_old);

                X = GM1\tmpF;
                Y = mxfem1d_F2L( obj.p, X);
                obj.E(:, j) = 2*Y-obj.E(:, j-1);

                indl = (nx_l)*(obj.p+1)+1:(nx_l)*(obj.p+1)+obj.p+1;
                indr = (nx_r)*(obj.p+1)+1:(nx_r)*(obj.p+1)+obj.p+1;
                yL(1, j) = exp(-1i*Xi1_l)*IMb_l*obj.E(indl,j);
                yR(1, j) = exp(-1i*Xi1_r)*IMb_r*obj.E(indr,j);
                
                SR_old = SR_new;
                SL_old = SL_new;
                
                whl0_old = whl0_new;
                whr0_old = whr0_new;

                u = obj.WP_data.(obj.soln_name)(obj.x, j, obj.soln_params);
                obj.E_ana(:, j) = mxfem1d_FLT( obj.N, obj.LGLdata.FM, u);

                obj.eTBC(j) = mxfem1d_Norm2( obj.p, obj.N,...
                                             obj.E(:,j)-obj.E_ana(:, j))/norm0;
            end
        end


        function obj = TBC_CQ_error(obj)
            dx_l    = obj.sw(1);
            dx_r    = obj.sw(2);

            % verify support of the initial data
            u = obj.WP_data.(obj.soln_name)(obj.x, 1, obj.soln_params);
            
            supp = LSE.find_support(obj.x, u, obj.stol);

            if ((obj.domain(1)>supp(1)-dx_l) || (obj.domain(2)<supp(2)+dx_r))

                fprintf( 'Support of the initial data: [%0.5f, %0.5f]\n',...
                          supp(1), supp(2));
                
                error('Computational domain is inadmissible!')
            end
            if~isempty(obj.phi)
                supp = LSE.find_support(obj.x1, obj.phi, obj.stol);
                if ((obj.domain(1)>supp(1)-dx_l) || (obj.domain(2)<supp(2)+dx_r))
                    fprintf( 'Support of the static potential: [%0.5f, %0.5f]\n',...
                              supp(1), supp(2));
                    error('Computational domain is inadmissible!')
                end
            end
            U = mxfem1d_FLT( obj.N, obj.LGLdata.FM, u);
            X = zeros(obj.nr_fem, 1);
            Y = zeros(obj.nr_lgl, 1);

            E = zeros(obj.nr_lgl, 1);
            E(:,1) = U;
            obj.eTBC = zeros(1, obj.Nt);
            norm0 = mxfem1d_Norm2( obj.p, obj.N, U);

            xl  = obj.domain(1) + dx_l - obj.tnu;
            xr  = obj.domain(2) - dx_r - obj.tnu;

            yL = [0];
            yR = [0];

            SL = 0;
            SR = 0;

            whl0_old = 0;
            whr0_old = 0;

            SR_old  = 0;
            SL_old  = 0;

            for j=2:obj.Nt,

                nu_tmp   = 0.5*(obj.tnu(j)   + obj.tnu(j-1));
                beta_tmp = 0.5*(obj.tbeta(j) + obj.tbeta(j-1));
                B_tmp    = 0.5*(obj.tB(j)    + obj.tB(j-1));

                whl = LSE.quadW_CQ(j, dx_l - obj.tnu(j), obj.dt, obj.tol);
                whr = LSE.quadW_CQ(j, dx_r + obj.tnu(j), obj.dt, obj.tol);

                tmpV = 0.5*(obj.V(j-1) + obj.V(j));
                GM1 = obj.GPMM*tmpV + obj.GM;

                xl_tmp  = 0.5*(xl(j) + xl(j-1));
                xr_tmp  = 0.5*(xr(j) + xr(j-1));

                [nx_l, xi_l] = LSE.fem_index(obj.N, obj.domain, xl_tmp);
                [nx_r, xi_r] = LSE.fem_index(obj.N, obj.domain, xr_tmp);

                IMb_l = mxLGLdataIM( obj.p, xi_l);
                IMb_r = mxLGLdataIM( obj.p, xi_r);
                 
                whl0_new = whl(1)*exp(1i*obj.tbeta(j)*( -dx_l + obj.tnu(j)));
                whr0_new = whr(1)*exp(1i*obj.tbeta(j)*(  dx_r + obj.tnu(j)));

                IMvb_l = 0.5*(whl0_new+whl0_old)*LSE.Lobatto(IMb_l);
                IMvb_r = 0.5*(whr0_new+whr0_old)*LSE.Lobatto(IMb_r);

                GM1(1, 1)      = GM1(1, 1)  + 1i*beta_tmp;
                GM1(1, nx_l+1) = GM1(1, nx_l+1) + IMvb_l(1);
                GM1(1, nx_l+2) = GM1(1, nx_l+2) + IMvb_l(2);
                
                shift = obj.N+1 + nx_l*(obj.p-2);

                for k=3:obj.p+1
                    GM1(1, shift+k) = GM1( 1, shift+k) + ...
                                           IMvb_l(k);
                end

                GM1(obj.N+1, obj.N+1) = GM1(obj.N+1, obj.N+1)  - 1i*beta_tmp;
                GM1(obj.N+1, nx_r+1)  = GM1(obj.N+1, nx_r+1) + IMvb_r(1);
                GM1(obj.N+1, nx_r+2)  = GM1(obj.N+1, nx_r+2) + IMvb_r(2);
                
                shift = obj.N + 1 + nx_r*(obj.p-2);
                
                for k=3:obj.p+1
                    GM1(obj.N+1, shift+k) = GM1( obj.N+1, shift+k) + ...
                                             IMvb_r(k);
                end

                Pvb = mxfem1d_PrjL2F( obj.p, E);
                F   = -1i*obj.r*obj.J*Pvb;


                whl_tmp = fliplr(whl((2:j)));
                whr_tmp = fliplr(whr((2:j)));

                Xi1_l = obj.tbeta(j)*xl(j)         - obj.tB(j); 
                Xi_l  = obj.tbeta(j)*obj.domain(1) - obj.tB(j); 

                Xi1_r = obj.tbeta(j)*xr(j)         - obj.tB(j); 
                Xi_r  = obj.tbeta(j)*obj.domain(2) - obj.tB(j); 

                SL_new = exp(1i*Xi_l)*sum(whl_tmp.*yL);
                SR_new = exp(1i*Xi_r)*sum(whr_tmp.*yR);


                tmpF = F;
                tmpF(1)       = F(1)       - 0.5*(SL_new + SL_old);
                tmpF(obj.N+1) = F(obj.N+1) - 0.5*(SR_new + SR_old);

                X = GM1\tmpF;
                Y = mxfem1d_F2L( obj.p, X);
                E = 2*Y - E;

                indl = (nx_l)*(obj.p+1)+1:(nx_l)*(obj.p+1)+obj.p+1;
                indr = (nx_r)*(obj.p+1)+1:(nx_r)*(obj.p+1)+obj.p+1;
                yL(1, j) = exp(-1i*Xi1_l)*IMb_l*E(indl,1);
                yR(1, j) = exp(-1i*Xi1_r)*IMb_r*E(indr,1);
                
                SR_old = SR_new;
                SL_old = SL_new;
                
                whl0_old = whl0_new;
                whr0_old = whr0_new;

                u = obj.WP_data.(obj.soln_name)(obj.x, j, obj.soln_params);
                U = mxfem1d_FLT( obj.N, obj.LGLdata.FM, u);
                obj.eTBC(j) = mxfem1d_Norm2( obj.p, obj.N, E-U)/norm0;
            end
        end
    end


    methods (Static=true, Access=private)

        function [tnu, tbeta, tB] = auxilary_func(dt, V, T_supp)
            % dnu=-2*beta, dbeta = V(t)
            % nu = -2partial^{-2}V(t)
            Nt = length(V);

            tbeta = zeros(1, Nt);
            tnu   = zeros(1, Nt);
            for i=2:Nt
                tbeta_old = tbeta(i-1);
                V_tmp = 0.5*(V(i)+V(i-1))*dt;
                tbeta(i)  = tbeta_old + V_tmp;
                tbeta_tmp = tbeta_old + 0.5*V_tmp;
                tnu(i)    = tnu(i-1) -2*dt*tbeta_tmp;
            end
            if(nargin>2)
                N_supp = ceil(T_supp/dt);
                if(Nt<N_supp)
                    error('Insuficient number of time-steps!');
                end
                tbeta_T = tbeta(N_supp);
                tbeta  = tbeta-tbeta_T;
                tnu    = tnu +2*tbeta_T*(0:Nt-1)*dt;
            end

            tB = zeros(1, Nt);
            for i=2:Nt
                tB(i)  = tB(i-1) + 0.5*dt*(tbeta(i)^2+tbeta(i-1)^2);
            end

        end


        % =========================================================
        % Convolution Quadrature, Trapezoidal rule
        % =========================================================
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
        % =========================================================

        function [nx, xi] = fem_index(N, domain, x )
            x_l = domain(1);
            x_r = domain(2);
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

    methods(Static)
        % =========================================================
        % Predefined potential functions
        % =========================================================
        function y = tpotential_CW(t, params)
            theta = params.mup*t;
            y     = params.g_0*cos(theta);
        end

        function y = tpotential_Gaussian(t, params)
            theta = params.mup*t;
            t0    = params.t_0;
            sigma = params.sigma;

            y = params.g_0*cos(theta).*exp(-(t-t0).^2/(sigma^2));
        end
        % =========================================================

        function supp = find_support(x, vec, tol)
        % find support of a function
            x_supp = x(abs(vec)>tol);
            supp = [min(x_supp), max(x_supp)];
        end
    end
end
