function superM = fem1d_superM( p, P, IMi, quadW)

sym i;
B = sym(zeros(p+1,p+1));
tmpB = sym([  1/2    1/2;...
             -1/2    1/2]);
B(1:2,1:2) = tmpB;
for i=3:p+1,
    B(i,i)      = 1/sqrt(2*(2*i-3));
    B(i-2,i)    = -1/sqrt(2*(2*i-3));
end

M = IMi*double(B);
superM = zeros(P+1, 3+(p-1)*(p+4)/2);

superM(:,1) = M(:,1).*M(:,1).*quadW;
superM(:,2) = M(:,1).*M(:,2).*quadW;
superM(:,3) = M(:,2).*M(:,2).*quadW;
for k=1:p-1,
    superM(:,k+3) = M(:,1).*M(:,k+2).*quadW;
    superM(:,k+3+p-1) = M(:,2).*M(:,k+2).*quadW;
end
shiftIn = 2*(p-1)+3;
for k=1:p-1,
    In = (1:k);
    ones_ = ones(1,k);
    superM( :, In+shiftIn) = M( :, In+2).*M(:,(k+2)*ones_).*quadW( :, ones(1,k));
    shiftIn = shiftIn+k;
end 
