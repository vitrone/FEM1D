function u = ILT(p, N, P, IM, v)
% N - number of (FEM) elements 
% p - highest degree of the Legendre polynomials involved
p1 = p+1;


u   =   zeros(P*N+1,1);
for i=0:N-1,
    tmpu            = IM*v(i*p1+1:i*p1+p1); 
    u(i*P+1:i*P+P)  = tmpu(1:end-1);
end
u(N*P+1) = tmpu(end);
