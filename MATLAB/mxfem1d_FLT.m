function mxfem1d_FLT
% mxfem1d_FLT Forward Legendre Transfrom
% U = mxfem1d_FLT( N, FM, u) computes forward legendre transfrom U of u. Here N is
% the number of finite elements, FM is the forward transform matrix.
%
% Example:
%       x_l = -5;
%       x_r =  5;
%       N   =  10;
%       p   =  4;
%       poly_fun = @(x) x.^p;
%       tol = 1e-9;
%       [xir, quadW, IM, FM] = mxLGLdataLT2( p, tol);
%       x = ref2mesh (xir, N, [x_l x_r]);
%       
%       u = poly_fun(x);      % a column vector
%       U = mxfem1d_FLT( N, FM, u);


