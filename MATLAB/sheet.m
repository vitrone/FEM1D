% symbolic calculations:
% =====================
clear
syms g_0 g_1 mu t beta ibeta nu g
g = g_0+g_1*cos(mu*t);
beta = int(g,'t');
nu = -2*int(beta,'t');
val_nu = subs(nu,{t},{0});
nu = nu - val_nu;
ibeta2 = int(beta^2,'t');
val_ibeta2 = subs(ibeta2,{t},{0});
ibeta2 = ibeta2-val_ibeta2;