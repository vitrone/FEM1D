%------------------------------------------------------------------------------%
% Name        : test_jacobi.m                                                  %
% Author      : Vishal Vaibhav                                                 %
%                                                                              %
% Description : Test numerical implementation of fractional operators based on %
% Jacobi polynomials using Pade approximation method.                          %
%                                                                              %
% History: Creation 23 Aug 2013                                                %
%                                                                              %
%------------------------------------------------------------------------------%

classdef test_jacobi < matlab.unittest.TestCase

    methods(Test)
        function test_partial_fraction_coeff(testCase)
        % compare the computed values with the analytically computed values 
        % for frac_pow=1/2;
            order       = 50; 
            frac_pow    = 1/2;
            K           = (1:order); 
            KK          = (0:order);
            theta       = K*pi/(2*order+1);
            D_ana       = ((tan(theta)).^2);
            beta0       = 2*order+1;
            beta        = -(2/(2*order+1))*(sin(theta)).^2./((cos(theta)).^4);
            N_ana       = ([beta0, beta]);
           
            % Numerically computed values
            tol = 1e-9;
            [N, D] = pf_pade_exp(order, frac_pow, tol);
            eN = norm(N_ana-N)/norm(N_ana);
            eD = norm(D_ana-D)/norm(D_ana);
            % fprintf('Numerical Error (eN, eD): (%0.12g, %0.12g )\n', eN, eD );

            testCase.verifyLessThan( [eN, eD], [tol, tol] );

        end
        %----------------------------------------------------------------------%
        % Following tests are based on results from 
        % Miller and Rose, Fractional Calculus
        %----------------------------------------------------------------------%

        function test_fractional_operator1(testCase)
        % function: t^(alpha), alpha >-1
        % D^(nu)t^(alpha) = [gamma(1+alpha)/gamma(1+alpha-nu)]t^(-nu+alpha), 
        % -1<nu<0
           fun = @(t, alpha) t.^alpha;
           alpha = 2.1;
           
           dt = 1e-4;
           Nt = 10000+1;
           r  = 2/dt;
           Mp = 50;
           nu = -0.6;
           
           tol = 1e-9;
           [Num Denom] = pf_pade_exp( Mp, nu, tol);
           b0    = Num(1);
           b_nu  = Num(2:end);
           G_nu  = r./(r+Denom);
           r_nu  = (b0+sum(b_nu.*G_nu)/r);
           bG_nu = 2*b_nu.*G_nu/r;

           psi  = 0;
           Psi  = 0;
           vphi = zeros(1, Mp);
           t    = dt/2;
           for i=1:Nt,
               psi  = fun(t, alpha);
               Psi  = r_nu*psi + sum(G_nu.*vphi);
               vphi = (2*G_nu-1).*vphi + bG_nu*psi;
               t    = t + dt;
           end
           K = (gamma(1+alpha)/gamma(1+alpha-nu));
           Psi_actual = K*fun( t-dt, alpha-nu);
           e = abs(Psi-Psi_actual)/abs(Psi_actual);

           testCase.verifyLessThan( e, tol );

        end
        
        function test_fractional_operator2(testCase)
        % function: t^(alpha), alpha >-1
        % D^(nu)t^(alpha) = [gamma(1+alpha)/gamma(1+alpha-nu)]t^(-nu+alpha), 
        % -1<nu<0
           fun = @(t, alpha) t.^alpha;
           alpha = 3.6;
           
           dt = 0.5e-4;
           Nt = 20000+1;
           r  = 2/dt;
           Mp = 50;
           nu = -0.6;
           
           tol = 1e-9;
           [Num Denom] = pf_pade_exp( Mp, nu, tol);
           b0    = Num(1);
           b_nu  = Num(2:end);
           G_nu  = r./(r+Denom);
           r_nu  = (b0+sum(b_nu.*G_nu)/r);
           bG_nu = 2*b_nu.*G_nu/r;

           psi  = 0;
           Psi  = 0;
           vphi = zeros(1, Mp);
           t    = dt/2;
           for i=1:Nt,
               psi  = fun(t, alpha);
               Psi  = r_nu*psi + sum(G_nu.*vphi);
               vphi = (2*G_nu-1).*vphi + bG_nu*psi;
               t    = t + dt;
           end
           K = (gamma(1+alpha)/gamma(1+alpha-nu));
           Psi_actual = K*fun( t-dt, alpha-nu);
           e = abs(Psi-Psi_actual)/abs(Psi_actual);

           testCase.verifyLessThan( e, tol );

        end
        
        function test_fractional_operator3(testCase)
        % function: e^(a*t)
        % D^(nu)e^(at) = a^(nu)*e^(a*t)*gammainc(at, -nu), -1<nu<0
        % gammainc - incomplete gamma function
           a = -1.2;
           dt = 1e-6;
           Nt = 1e6+1;
           r  = 2/dt;
           Mp = 50;
           nu = -0.6;
           
           tol = 1e-9;
           [Num Denom] = pf_pade_exp( Mp, nu, tol);
           b0    = Num(1);
           b_nu  = Num(2:end);
           G_nu  = r./(r+Denom);
           r_nu  = (b0+sum(b_nu.*G_nu)/r);
           bG_nu = 2*b_nu.*G_nu/r;

           psi  = 0;
           Psi  = 0;
           vphi = zeros(1, Mp);
           t    = dt/2;
           for i=1:Nt,
               psi  = exp(a*t);
               Psi  = r_nu*psi + sum(G_nu.*vphi);
               vphi = (2*G_nu-1).*vphi + bG_nu*psi;
               t    = t + dt;
           end
           % Psi = Psi + (t-dt)^(-nu)/gamma(1-nu);
           K = a^(nu)*exp(a*t);
           Psi_actual = K*gammainc(a*(t-dt), -nu );
           e = abs(Psi-Psi_actual)/abs(Psi_actual);

           testCase.verifyLessThan( e, 1e-5 );

        end

        function test_fractional_operator4(testCase)
        % function: cos(a*t)
        % D^(nu)cos(a*t) = t^(-nu)*h(1,[1-nu, 2-nu]/2; -(a*t)^2/4)/gamma(1-nu)],
        % -1<nu<0, where h refers to the hyper-geometric function
           a  = pi/2;
           dt = 1e-3;
           Nt = 1e3+1;
           r  = 2/dt;
           Mp = 50;
           nu = -0.6;
           
           tol = 1e-9;
           [Num Denom] = pf_pade_exp( Mp, nu, tol);
           b0    = Num(1);
           b_nu  = Num(2:end);
           G_nu  = r./(r+Denom);
           r_nu  = (b0+sum(b_nu.*G_nu)/r);
           bG_nu = 2*b_nu.*G_nu/r;

           psi  = 0;
           Psi  = 0;
           vphi = zeros(1, Mp);
           t    = dt/2;
           for i=1:Nt,
               psi  = cos(a*t);
               Psi  = r_nu*psi + sum(G_nu.*vphi);
               vphi = (2*G_nu-1).*vphi + bG_nu*psi;
               t    = t + dt;
           end
           K = (t-dt)^(-nu)/gamma(1-nu);
           Psi_actual = K*hypergeom(1, [(1-nu)/2, (2-nu)/2], -0.25*a^2*(t-dt)^2);
           e = abs(Psi-Psi_actual)/abs(Psi_actual);

           testCase.verifyLessThan( e, 1e-6 );

        end

    end
end
