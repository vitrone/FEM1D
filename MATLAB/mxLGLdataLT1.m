function mxLGLdataLT1
% mxLGLdataLT1 Data for Discrete Legendre Transforms
%     [zeros, quadW] = LGLdataLT1(degree, tolerance) 
%     Calculates the zeros of the LGL polynomial and the associated Gauss 
%     quadrature weights.
%     degree    - highest degree of the polynomials used
%     quadW     - Gauss quadrature weights corresponding to the LGL points
%     
%     Example: 
% 
%     f = @(x) exp(-16*x.^2);
%     [zeros quadW] = LGLdataLT1(60,1e-9);
%     F = f(zeros);   % a column vector
%     Integral = sum(F.*quadW) % integral on [-1,1], exact value = sqrt(pi/16)
