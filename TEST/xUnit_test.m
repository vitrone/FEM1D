clc
clear
set_pathvariables;
tc = xUnit_fem1d

%run(tc, 'test_interpolation')
%run(tc, 'test_interpolation2')
%run(tc, 'test_norm')
%run(tc, 'test2_norm')
run(tc, 'test_GP_ABC1a_CQ')
%run(tc, 'test_fem1d_L2F')
%run(tc, 'test_fem1d_sparse_GMM')

%x_l = -5;
%x_r =  5;
%N   =  200;
%p   =  4;);
%nr_LGL = 2*p+1;
%P =  nr_LGL-1;);
%
%fun = @(x) exp(-x.^2/2);
%pot = @(x) 1+0*x*1i;
%tol = 1e-9;
%
%[xi, quadW] = mxLGLdataLT1( P, tol);
%
%FM = mxLGLdataFM(p, xi);
%IM = mxLGLdataIM(p, xi);
%Q  = mxfem1d_quadM(quadW, IM);
%x  = ref2mesh (xi, N, [x_l x_r]);
%
%% u -> a column vector
%u = fun(x);      
%U = mxfem1d_FLT( N, FM, u);
%vb = mxfem1d_L2F(p, U);
%U2 = mxfem1d_F2L(p, vb);
%
%
%phi = pot(x);
%M = mxfem1d_sparse_GMM(p, Q, phi);
%
%u1 = u.*phi;
%U1 = mxfem1d_FLT( N, FM, u1);
%Pvb = mxfem1d_PrjL2F(p, U1);
%
%sM = M + transpose(M) - diag(diag(M));
%Pvb1 = sM*vb;
%
%err = (Pvb-Pvb1);
%
%e_relative = norm(err)/norm(Pvb)




