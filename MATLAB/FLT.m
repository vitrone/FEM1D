function v = FLT(p, N, P, FM, u)
% N - number of (FEM) elements 
% p - highest degree of the Legendre polynomial involved
% P - P+1 is the number of LGL points used for sampling
% Note that P>=p. If P>p then he transformation is truncated 
% after the polynomial degree p, the ouput vector is of length p+1 
% per FEM-element instead of P+1.  

p1 = p+1;

v = zeros(N*p1,1);
for i=0:N-1,
    tmpv = FM*u(i*P+1:i*P+P+1); 
    v(i*p1+1:i*p1+p1) = tmpv;
end
