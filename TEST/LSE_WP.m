classdef LSE_WP

% Analytical solutions of the Schrodinger equation
% iu_t + u_xx + xV(t)u = 0.
% derived from the known solution of the free Schrodinger equation
% iu_t + u_xx = 0.
    
    properties
        dt   % time-step
        Nt   % total number of points in t
        t
        % dnu=-2*beta, dbeta = V(t)
        % nu = -2partial^{-2}V(t)
        beta 
        nu
        B
        V
        params
    end


    methods (Static)

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

        function [nu, beta, B] = auxilary_func(dt, V)
            % dnu=-2*beta, dbeta = V(t)
            % nu = -2partial^{-2}V(t)
            Nt = length(V);

            beta = zeros(1, Nt);
            nu   = zeros(1, Nt);
            
            for i=2:Nt
                beta_old = beta(i-1);
                V_tmp = 0.5*(V(i)+V(i-1))*dt;
                beta(i)  = beta_old + V_tmp;
                beta_tmp = beta_old + 0.5*V_tmp;
                nu(i)    = nu(i-1) -2*dt*beta_tmp;
            end

            B = zeros(1, Nt);
            for i=2:Nt
                B(i) = B(i-1) + 0.5*dt*(beta(i)^2+beta(i-1)^2);
            end

        end

        function [coeff, A] = Hermite_poly(n)
            switch(n)
                case{0}
                    coeff = [1];
                    A = 1;
                case{1}
                    coeff = [2, 0];
                    A = 1/sqrt(2*sqrt(pi));
                case{2}
                    coeff = [4, 0, -2];
                    A = 1/sqrt(2*sqrt(pi));
                case{3}
                    coeff = [8, 0, -12, 0];
                    A = 1/sqrt(48*sqrt(pi));
                case{4}
                    coeff = [16, 0, -48, 0, 12];
                    A = 1/sqrt(384*sqrt(pi));
                otherwise
                    error('Higher degree polynomials are not available!')
            end
        end
    end

    methods
        function obj = LSE_WP(dt, Nt, params, V)
            % defaults
            obj.dt = dt;
            obj.Nt = Nt;
            obj.t  = (0:Nt-1)*dt;
            obj.V  = zeros(1, Nt);
            obj.params = params;

            if(nargin>3)
                obj.V = V;
            else
                switch(params.potential_name)
                    case{'CW'}
                        obj.V = LSE_WP.tpotential_CW(obj.t, params);
                    case{'GAUSSIAN'}
                        obj.V = LSE_WP.tpotential_Gaussian(obj.t, params);
                    otherwise
                        error('potential is unknown!')
                end
            end
            [obj.nu, obj.beta, obj.B] = LSE_WP.auxilary_func(dt, obj.V);
        end

        function y = Gaussian_WP(obj, x, j, sol_params)
            t_tmp = obj.t(j);
            c = sol_params.c;

            x1 = x + obj.nu(j) - c*t_tmp;
            Xi = 0.5*c*x1 + 0.25*c^2*t_tmp + obj.beta(j)*x - obj.B(j);

            C  = (1/(4*sol_params.a*t_tmp*1i+1));
            y  = (sol_params.A0)*sqrt(C)*exp(-sol_params.a*C*(x1).^2+1i*Xi);
        end

        function y = Gaussian_WP_mesh(obj, xx, JJ, sol_params)
            [x,J] = ndgrid( xx, JJ);
            t_tmp = obj.t(J);
            c = sol_params.c;

            x1 = x + obj.nu(J) - c*t_tmp;
            Xi = 0.5*c*x1 + 0.25*c^2*t_tmp + obj.beta(J).*x - obj.B(J);

            C  = (1./(4*sol_params.a*t_tmp*1i+1));
            y  = (sol_params.A0)*sqrt(C).*exp(-sol_params.a*C.*(x1).^2+1i*Xi);
        end


        function y = GaussHermite_WP(obj, x, j, sol_params)
            t_tmp = obj.t(j);
            c = sol_params.c;
            a = sol_params.a;

            x1 = x + obj.nu(j) - c*t_tmp;
            Xi = 0.5*c*x1 + 0.25*c^2*t_tmp + obj.beta(j)*x - obj.B(j);

            C  = (1/(4*a*t_tmp*1i+1));
            y_tmp  = (sol_params.A0)*sqrt(C)*exp(-a*C*(x1).^2+1i*Xi);
            [coeff, A] = LSE_WP.Hermite_poly(sol_params.degree);

            x1_tmp = sqrt(2*a)*x1/sqrt(1+(4*a*t_tmp)^2);
            h_tmp  = A*polyval(coeff, x1_tmp)/sqrt(2*a);
            e_tmp  = exp(-1i*sol_params.degree*atan(4*t_tmp*a));

            y = y_tmp.*e_tmp.*h_tmp;
        end
        
        function y = GaussHermite_WP_mesh(obj, xx, JJ, sol_params)
            [x,J] = ndgrid( xx, JJ);
            t_tmp = obj.t(J);
            c = sol_params.c;
            a = sol_params.a;

            x1 = x + obj.nu(J) - c*t_tmp;
            Xi = 0.5*c*x1 + 0.25*c^2*t_tmp + obj.beta(J).*x - obj.B(J);

            C  = (1./(4*a*t_tmp*1i+1));
            y_tmp  = (sol_params.A0)*sqrt(C).*exp(-a*C.*(x1).^2+1i*Xi);
            
            [coeff, A] = LSE_WP.Hermite_poly(sol_params.degree);
            disp(coeff);

            x1_tmp = sqrt(2*a)*x1./sqrt(1+(4*a*t_tmp).^2);
            h_tmp  = A*polyval(coeff, x1_tmp)/sqrt(2*a);
            e_tmp  = exp(-1i*sol_params.degree*atan(4*t_tmp*a));

            y = y_tmp.*e_tmp.*h_tmp;
        end

    end
end
