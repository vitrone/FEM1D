function x = ref2mesh(xi,N,W)
x_l = W(1);
x_r = W(2);
dx  = (x_r-x_l)/N;
J   =   dx/2;
% ref_map = @(xi,J,xm) J*xi+xm;
p   =   length(xi)-1;
x   =   zeros(p*N+1,1);
for i=0:N-1,
    %     xm=x_l+J
    tmpx                = J*xi+x_l+(2*i+1)*J;
    x(i*p+1:i*p+p,1)    = tmpx(1:end-1);
end
x(p*N+1) = x_r;
