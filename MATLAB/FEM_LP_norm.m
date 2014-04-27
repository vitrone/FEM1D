function norm = FEM_LP_norm( p, N, u)
norm = 0;
K     =   transpose(0:p);
for i=0:N-1, 
    norm = norm...
        +sum(abs(u(i*(p+1)+1:i*(p+1)+p+1)).^2.*(2./(2*K+1)));
end
norm = sqrt(norm);