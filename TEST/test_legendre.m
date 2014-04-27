classdef test_legendre < matlab.unittest.TestCase
% test = test_legendre
% run(test)

    methods(Test)
        function test_quadrature(testCase)
            p = 5;
            poly_fun = @(x) x.^(2*p+1);
            tol = 1e-9;
            [xir, quadW] = LGLdataLT1(p, tol);
            F = poly_fun(xir);      % a column vector
            integral = sum(F.*quadW);
            % fprintf('Value of the integral = %0.12g\n', integral);
            testCase.verifyLessThan( integral, tol );

        end
        function test_quadrature2(testCase)
            p = 50;
            fun = @(x) exp(-16*x.^2);
            tol = 1e-9;
            [xir, quadW] = LGLdataLT1(p, tol);
            F = fun(xir);      % a column vector
            integral = sum(F.*quadW);
            expected_val = sqrt(pi/16); 
            e = abs(integral-expected_val)/abs(expected_val);
            % fprintf('Value of the integral = %0.12g\n', integral);
            res_tol = 1e-6;
            testCase.verifyLessThan( e, res_tol );

        end
        function test_interpolation(testCase)
            % Legendre transform is exact for polynomials
            p = 5;
            poly_fun = @(x) x.^p;
            tol = 1e-9;
            [xir, quadW, IM, FM] = LGLdataLT2(2*p, tol);
            F = poly_fun(xir);      % a column vector
            U = FM*F;               % discrete Legendre transform

            xi = (-1:0.001:1)';
            P = length(xi)-1;
            IM1 = LGLdataIM( P, 2*p, xi);
            X = IM1*U;              % inverse discrete Legendre transform
            F = poly_fun(xi);
            e = norm(F-X)/norm(F);
            testCase.verifyLessThan( e, tol );

        end
        function test_anti_aliasing(testCase)
            % Test anti-aliasing 
            p = 60;
            d = 3;
            tol = 1e-9;
            fun = @(x) exp(-16*x.^2);
            [xir, quadW, IM, FM] = LGLdataLT2(p, tol);
            F = fun(xir);           % a column vector
            U = FM*F;               % discrete Legendre transform

            xi = (-1:0.001:1)';
            P = length(xi)-1;
            IM0 = LGLdataIM( P, p, xi);
            F2  = (fun(xi)).^d;

            % FIRST: Find the error without anti-aliasing
            P = p;
            [xir, quadW] = LGLdataLT1(P, tol);
            IM = LGLdataIM( P, p, xir);
            X  = IM*U;               % inverse discrete Legendre transform
            X2 = X.^d; % Note that the basis would require maximum degree=0.5*d*p
            FM = LGLdataFM(p, P, xir);
            U2 = FM*X2;             % discrete Legendre transform

            Xtmp = IM0*U2;          % inverse discrete Legendre transform
            e0   = norm(F2-Xtmp)/norm(F2);

            % Second: Find the error with anti-aliasing
            P = ceil(d*p*0.5)+1;
            [xir, quadW] = LGLdataLT1(P, tol);
            IM = LGLdataIM( P, p, xir);
            X  = IM*U;               % inverse discrete Legendre transform
            X2 = X.^d; % Note that the basis would require maximum degree=0.5*d*p
            FM = LGLdataFM(p, P, xir);
            U2 = FM*X2;             % discrete Legendre transform

            Xtmp = IM0*U2;          % inverse discrete Legendre transform
            e    = norm(F2-Xtmp)/norm(F2);
            % Finally compare the two errros!
            testCase.verifyLessThan( e, e0 );
        end
    end
end
