% x = ref2mesh( xr, N, [x_l x_r]); m
% [xr quadW IM FM] = LGLdataLT2( p, tol); mex
% [xr quadW] = LGLdataLT1(p, tol); mex
% FM = LGLdataFM( p, P, xr); mex
% IM = LGLdataIM( P, p, xr); mex

% w = FLT(p, N, P, FM, u); m
% u = ILT( p, N, P, IM, w); m
% w = fem_flt( p, N, P, FM, u); mex
% u = fem_ilt( p, N, P, IM, w); mex

% [Pv Pb] = prjlp2fem( p, N, w); mex
% [Pv Pb] = prjLP2FEM_ShapeFunc( p, N, w); m

% [v b] = LP2FEM_ShapeFunc( p, N, w); m
% w = FEM_ShapeFunc2LP( p, N, vv, bb); m

% [v b] = lp2fem( p, N, w); mex 
% w = fem2lp( p, N, v, b); m
% norm = FEM_LP_norm2( p, N, u); m
% norm = femlp_norm2( p, N, u); mex


% [N,D] = pf_pade_exp( order, frac_pow, tol); mex

% superM = fem1d_superM( p, P, IMi, quadW); m
% [row col nnz] = fem1d_crgmmind( p, N); m
% [gpmm lgpmm] = fem1d_crgmm( p, N, P, pot, superM, nnz); m
% [gpmm lgpmm] = fem1d_gmm( p, N, P, pot, superM, nnz); mex

% [row col gsm nnz] = fem1d_crgsm( p, N); m
