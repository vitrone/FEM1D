function mxfem1d_sparse_GMM
% M = fem1d_sparse_GMM(p, Q, phi) computes the global mass matrix for the given
% potential vector phi with the help of a precomputed quadrature matrix Q.
% p: maximum degree of the legendre poly in each FEM-elements
%
% Example: 
% x_l = -5;
% x_r =  5;
% N   =  20;
% p   =  4;
% nr_LGL = 2*p+1;
% P =  nr_LGL-1;
% 
% fun = @(x) exp(-x.^2/2);
% pot = @(x) 1+x.^2;
% tol = 1e-9;
% 
% [xi, quadW] = mxLGLdataLT1( P, tol);
% 
% FM = mxLGLdataFM(p, xi);
% IM = mxLGLdataIM(p, xi);
% Q  = mxfem1d_quadM(quadW, IM);
% x  = ref2mesh (xi, N, [x_l x_r]);
% 
% % u -> a column vector
% u = fun(x);      
% U = mxfem1d_FLT( N, FM, u);
% 
% phi = pot(x);
% M = mxfem1d_sparse_GMM(p, Q, phi);
