function pf_pade_exp
 
% Finds the Mth order diagonal Pade approximation of exponential functions 
% of the form s^a,|a|<1 (frac_pow=a) in the form of partial fractions. 
% Computes zeros of J^{(-a,a)}_M to find the numerator and the denominator 
% coefficients: 
% s^a=\beta_0+\sum_{k=0}^M\frac{\beta_k}{s+\gamma_k}.
% 
% [ N, D] = pf_pade_exp( order, frac_pow, tol);
% order - integer, order of the diagonal pade approximants
% frac_pow - Real number with modulus less than 1
% tol - tolerance
% N - numerator coefficients, order+1-by-1 array
% D - denomenator coefficients, order-by-1 array
% 
% 
% Example1:
% compare the computed values with the analytically computed values 
% for frac_pow=1/2;
% order = 50; frac_pow = 1/2;
% K = (1:order); KK = (0:order);
% theta = K*pi/(2*order+1);
% D_ana = (tan(theta)).^2;
% beta0 = 2*order+1;
% beta  = -(2/(2*order+1))*(sin(theta)).^2./((cos(theta)).^4);
% N_ana = transpose([beta0, beta]);
% 
% [N,D] = pf_pade_exp(order, frac_pow, 1e-9);
% 
% subplot(2,1,1), semilogy(KK, abs(N), 'r.', KK, abs(N_ana),'bo');
% title('Numerator Coefficients')
% legend('analytical', 'numerical','Location', 'EastOutside');
% subplot(2,1,2), semilogy(K, D, 'r.', K, D_ana,'bo')
% title('Denomenator Coefficients')
% legend('analytical', 'numerical','Location', 'EastOutside');
