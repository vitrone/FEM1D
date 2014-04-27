function [row col nnz] = fem1d_crgmmind( p, N)
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
%   row - contains the row indices
%   col - contains the column indices
%   gpmm - contains values of the nonzero elements as described by row 
%         and col vectors

nnz = N*p*(p+3)/2+1;

row = zeros( nnz, 1);
col = zeros( nnz, 1);

s = 1;
col(s) = 1; row(s) = 1;
s = 2;
col(s) = 2; row(s) = 1;
s = 3;
col(s) = 2; row(s) = 2;
for k=3:N,
    s = s+1;
    col(s) = k; row(s) = k-1;
    s = s+1;
    col(s) = k; row(s) = k; 
end
k = N+1;
s = s+1;
col(s) = k; row(s) = k-1;

s = s+1;
col(s) = k; row(s) = k;
st1 = N+1;
for i=1:N,
    for k = 1:p-1,
        s = s+1;
        col(s) = st1+k; row(s) = i;
        s = s+1;
        col(s) = st1+k; row(s) = i+1;
        for l=1:k,
            s = s+1;
            col(s) = st1+k; row(s) = l+st1;
        end
    end
    st1 = st1+p-1;
end