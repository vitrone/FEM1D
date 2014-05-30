function mxLGLdataLT2
% LGLdataLT2 Data for Discrete Legendre Transforms
%     [zeros, quadW, IM, FM] = LGLdataLT(degree, tolerance) calculates neccessary 
%     data for discrete Legendre transform on Legendre-Gauss-Lobbato (LGL) 
%     points. degree is the highest degree of the Legenre polynomial involved. 
%     tolerance is need for the root finding algorithm to compute the degree+1 
%     LGL points. The output arguments are:
%     zeros  - degree+1 LGL points
%     quadW  - quadrature weights corresponding to the LGL points
%     IM and FM - backward and forward transform matrices, respectively
%     
%     Example: 
%     f = @(x) exp(-16*x.^2);
%     [zeros, quadW, IM, FM] = LGLdataLT2(60,1e-9);
%     F = f(zeros);   % a column vector
%     U = FM*F;       % discrete Legendre transform
%     X = IM*U;       % inverse discrete Legendre transform
%     error = abs(F-X);
%     semilogy(zeros,error)
%     I = sum(F.*quadW) % integral on [-1,1], exact value = sqrt(pi/16)
