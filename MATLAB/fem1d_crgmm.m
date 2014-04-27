function [gpmm lgpmm] = fem1d_crgmm( p, N, P, pot, superM, nnz)
% This function performs Gauss quadrature on LGL points to determine the 
% upper triangular half of the symmteric mass matrix for a given potential 
% function.
% 
% p - highest degree of the Legendre polynomial involved
% P+1 - number of LGL points used for integration 
% M - (P+1)-by-(p+1) matrix containing values of the basis function
% on LGL points along with the quadrature weights.
% pot - vector containing values of the potential function on LGL points
% 
% nnz - number of nonzero elements
% The following are column vectors:
%   gpmm - contains values of the nonzero elements as described by row 
%         and col indices vectors

Ni = P+1;
% nnz = N*p*(p+3)/2+1;

% row = zeros( nnz, 1);
% col = zeros( nnz, 1);
gpmm = zeros( nnz, 1);

K = (1:Ni);
st = (Ni-1);

s = 1;
% col(s) = 1; row(s) = 1;
gpmm(s) = sum(pot(K).*superM(:,1));
s = 2;
% col(s) = 2; row(s) = 1;
gpmm(s) = sum(pot(K).*superM(:,2));
s = 3;
% col(s) = 2; row(s) = 2;
gpmm(s) = sum(pot(K).*superM(:,3))+sum(pot(st+K).*superM(:,1));
tmpK = K;
for k=3:N,
    s = s+1;
%     col(s) = k; row(s) = k-1;
    tmpK = st+tmpK;
    gpmm(s) = sum(pot(tmpK).*superM(:,2));
    s = s+1;
%     col(s) = k; row(s) = k; 
    gpmm(s) = sum(pot(tmpK).*superM(:,3))+sum(pot(st+tmpK).*superM(:,1));
    
end
k = N+1;
s = s+1;
% col(s) = k; row(s) = k-1;
tmpK = st+tmpK;
gpmm(s) = sum(pot(tmpK).*superM(:,2));

s = s+1;
% col(s) = k; row(s) = k;
gpmm(s) = sum(pot(tmpK).*superM(:,3));

tmpK = K;
st1 = N+1;
for i=1:N,
    for k = 1:p-1,
        s = s+1;
%         col(s) = st1+k; row(s) = i;
        gpmm(s) = sum(pot(tmpK).*superM(:,3+k));
        s = s+1;
%         col(s) = st1+k; row(s) = i+1;
        gpmm(s) = sum(pot(tmpK).*superM(:,3+k+p-1));
        tmp = 3+2*(p-1)+k*(k-1)/2;
        for l=1:k,
            s = s+1;
%             col(s) = st1+k; row(s) = l+st1;
            gpmm(s) = sum(pot(tmpK).*superM(:,tmp+l));
        end
    end
    st1 = st1+p-1;
    tmpK = st+tmpK;
end


lgpmm = gpmm;
s = 1;
lgpmm(s) = 0;
for k=2:N+1,
    s = s+2;
    lgpmm(s) = 0; 
end
for i=1:N,
    for k = 1:p-1,
        s = s+2+k;
        lgpmm(s) = 0;
    end
end





