function [row col gsm nnz] = fem1d_crgsm( p, N)
% In future this function will compute the stiffness matrix for general
% Sturm-Liouville problems

nnz = N*(p+2)+1;

row = zeros( nnz, 1);
col = zeros( nnz, 1);
gsm = zeros( nnz, 1);

s = 1;
col(s) = 1; row(s) = 1;
gsm(s) = 0.5;
s = 2;
col(s) = 1; row(s) = 2;
gsm(s) = -0.5;
s = 3;
col(s) = 2; row(s) = 1;
gsm(s) = -0.5;
s = 4;
col(s) = 2; row(s) = 2;
gsm(s) = 1;
s = 5;
col(s) = 2; row(s) = 3;
gsm(s) = -0.5;
for k=3:N,
    s = s+1;
    col(s) = k; row(s) = k-1;
    gsm(s) = -0.5;
    s = s+1;
    col(s) = k; row(s) = k; 
    gsm(s) = 1;
    s = s+1;
    col(s) = k; row(s) = k+1; 
    gsm(s) = -0.5;
end
k = N+1;
s = s+1;
col(s) = k; row(s) = k-1;
gsm(s) = -0.5;

s = s+1;
col(s) = k; row(s) = k;
gsm(s) = 0.5;

st = N+1;
for i=1:N,
    st1 = (i-1)*(p-1);
    for k = 1:p-1,
        s = s+1;
        col(s) = st+st1+k; row(s) = col(s);
        gsm(s) = 1;
    end
end