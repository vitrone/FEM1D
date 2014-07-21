clc
clear
%%
I = sqrt(-1);
tol = 1e-9; % tolerance for LGL points
N   = 400;  % Number of finite elements
p   = 4;    % highest degree of polynomials
[xr, quadW, IM, FM] = mxLGLdataLT2( p, tol);
% xr: LGL points
% FM: matrix for forward transform
% IM: matrix for inverse transform

P = 2*p; % For the purpose of numerical integration
% Degree of the Legendre polynomials involved --> P
[xri, quadWi] = mxLGLdataLT1( P, tol);
FMi = mxLGLdataFM( p, xri);
IMi = mxLGLdataIM( p, xri);

dt  =  0.5e-3;
r   =  2/dt;
Nt  =  2^14;
t   = (0:Nt-1)*dt;
x_r =  10;
x_l = -10;

h = (x_r-x_l)/N;
dx_l = 1.03;
dx_r = 1.03;

x   =  ref2mesh(xr,N,[x_l x_r]);
A0 = 1;
a  = 1/2; 
cl  = -5;
cr  = 5;

g0 = 1;
mup = 2*pi;

profile = @(x,t)LSE.Gaussian_WP3(x,t, A0, a, cr, g0, mup);
% profile = @(x,t)LSE.Gaussian_WP3(x,t, A0, a, cr, g0, mup)+...
%                 LSE.Gaussian_WP3(x,t, A0, a, cl, g0, mup);

nu   = @(t)LSE.nu(t, g0, mup);
beta = @(t)LSE.beta(t, g0, mup);
B    = @(t)LSE.B(t, g0, mup);
V    = @(t)LSE.V(t, g0, mup);

Data = struct( 'N'     , N,...
               'p'     , p,...
               'P'     , P,...
               'FM'    , FM,...
               'IM'    , IM,...
               'FMi'   , FMi,...
               'IMi'   , IMi,...
               'quadWi', quadWi,...
               'xri'   , xri,...
               'dt'    , dt,...
               'Nt'    , Nt,...
               'x_r'   , x_r,...
               'x_l'   , x_l,...
               'dx_l'  , dx_l,...
               'dx_r'  , dx_r,...
               'x'     , x,...
               'hfunc' , profile,...
               'nu'    , nu,...
               'B'     , B,...
               'beta'  , beta,...
               'V'     , V,...
               'tol'   , tol );

%%
% eTBC = LSE.TBC_CQ_error3(Data);
% logerr_plot(t, eTBC)
%%
E = LSE.TBC_CQ_evol3(Data);
%%
xi  = (-1:0.1:1);
x1  = ref2mesh( xi, N, [x_l x_r]);
tmpIM = mxLGLdataIM( p, xi);

% K = (1:20:Nt);
PSI = mxfem1d_ILT( N, tmpIM, E);

%%
options =   struct('levels',(-8:1),...
                'magnification_factor',1e0,...
                'tolerance', tol);
% contour_plot(x1,t,PSI,options)
K1 = (1:500:Nt);
curve_3d(x1,t,PSI,K1);
% data to be plotted --> x,t,PSI
